import math
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from calibration.refresh_saved_calibration import ee_to_ddg, fit, scores


class CalibrationRefreshTests(unittest.TestCase):
    def test_ee_conversion_scales_with_temperature_and_preserves_sign(self):
        old = ee_to_ddg(87, 298.15)
        self.assertAlmostEqual(ee_to_ddg(87, 333) / old, 333 / 298.15)
        self.assertAlmostEqual(ee_to_ddg(-87, 333), -ee_to_ddg(87, 333))
        for ee, temperature in [(100, 333), (-100, 333), (87, 0), (math.nan, 333)]:
            with self.assertRaises(ValueError):
                ee_to_ddg(ee, temperature)

    def test_refit_scales_coefficients_and_errors_but_preserves_r2(self):
        train_x, train_y = [1, 2, 3], [1.1, 2.4, 2.8]
        test_x, test_y = [1.5, 2.5], [1.4, 2.7]
        slope, intercept = fit(train_x, train_y)
        original = scores(test_y, [slope * x + intercept for x in test_x])
        factor = 333 / 298.15
        new_slope, new_intercept = fit(train_x, [v * factor for v in train_y])
        updated = scores([v * factor for v in test_y],
                         [new_slope * x + new_intercept for x in test_x])
        self.assertAlmostEqual(new_slope, factor * slope)
        self.assertAlmostEqual(new_intercept, factor * intercept)
        self.assertAlmostEqual(updated['MAE'], factor * original['MAE'])
        self.assertAlmostEqual(updated['R2'], original['R2'])

    def test_test_mean_null_has_zero_predictive_r2(self):
        metric = scores([1, 2, 3], [2, 2, 2])
        self.assertEqual(metric['R2'], 0)
        self.assertIsNone(metric['rho'])

    def test_family_commands_reproduce_complete_outputs_from_another_directory(self):
        source = Path(__file__).resolve().parents[1] / 'calibration'
        for family, folder in [('ni', 'nickel_epoxide_arylation'), ('pd', 'palladium_enone_addition')]:
            with self.subTest(family=family), tempfile.TemporaryDirectory() as tmp:
                base = Path(tmp) / 'calibration'
                target = base / folder
                target.mkdir(parents=True)
                shutil.copy2(source / 'refresh_saved_calibration.py', base)
                shutil.copy2(source / folder / 'calibrate.py', target)
                for path in (source / folder).glob('*.csv'):
                    shutil.copy2(path, target)
                expected = {p.name: p.read_bytes() for p in target.glob('*.csv')}
                for suffix in ('calibration_plot_data', 'calibration_plot_points', 'held_out_predictions'):
                    (target / f'{family}_{suffix}.csv').unlink()
                # A stale summary must be rebuilt with the null-model columns intact.
                summary = target / 'calibration_summary.csv'
                with summary.open(newline='') as handle:
                    row = next(csv.DictReader(handle))
                row['training_mean_null_test_MAE_mean'] = '-999'
                with summary.open('w', newline='') as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(row))
                    writer.writeheader()
                    writer.writerow(row)
                completed = subprocess.run([sys.executable, str(target / 'calibrate.py')],
                                           cwd=tmp, capture_output=True, text=True)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                for name, content in expected.items():
                    self.assertEqual((target / name).read_bytes(), content, name)
                metadata = json.loads((target / 'calibration_refit.json').read_text())
                self.assertEqual(metadata['n_splits'], 100)
                self.assertEqual(metadata['seed'], 20260810)
                self.assertTrue(metadata['temperature_source'].startswith('Reaction study, DOI'))


if __name__ == '__main__':
    unittest.main()
