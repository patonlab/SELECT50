import argparse
import importlib.util
import os
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / 'examples/calculations/opt_mace_omol.py'
spec = importlib.util.spec_from_file_location('calculation_example', SCRIPT)
example = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example)


class TrajectorySelectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.input_dir = self.root / 'input'
        self.run_dir = self.root / 'run'
        self.input_dir.mkdir()
        self.run_dir.mkdir()
        self.args = argparse.Namespace(xyz=self.input_dir / 'chosen.xyz')

    def trajectory(self, directory, name, timestamp):
        path = directory / name
        path.write_text('trajectory fixture')
        os.utime(path, (timestamp, timestamp))
        return path

    def test_expected_run_trajectory_wins_over_newer_unrelated_files(self):
        expected = self.trajectory(self.run_dir, 'chosen_optim.xyz', 1000)
        self.trajectory(self.run_dir, 'other_optim.xyz', 2000)
        self.trajectory(self.input_dir, 'chosen_optim.xyz', 3000)
        self.assertEqual(example.find_geometric_trajectory(self.args, self.run_dir), expected)

    def test_expected_input_trajectory_wins_over_newer_fallback(self):
        expected = self.trajectory(self.input_dir, 'chosen_optim.xyz', 1000)
        self.trajectory(self.run_dir, 'other_optim.xyz', 2000)
        self.assertEqual(example.find_geometric_trajectory(self.args, self.run_dir), expected)

    def test_fallback_requires_one_unique_trajectory(self):
        self.assertIsNone(example.find_geometric_trajectory(self.args, self.run_dir))
        fallback = self.trajectory(self.run_dir, 'legacy_optim.xyz', 1000)
        self.assertEqual(example.find_geometric_trajectory(self.args, self.run_dir), fallback)
        self.trajectory(self.input_dir, 'other_optim.xyz', 2000)
        with self.assertRaisesRegex(ValueError, 'Ambiguous'):
            example.find_geometric_trajectory(self.args, self.run_dir)


if __name__ == '__main__':
    unittest.main()
