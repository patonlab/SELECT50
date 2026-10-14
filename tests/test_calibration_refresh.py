import math
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


if __name__ == '__main__':
    unittest.main()
