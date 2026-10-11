import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

try:
    import numpy as np
    from scripts import calculate_heavy_atom_rmsd as rmsd
except ModuleNotFoundError as exc:
    if exc.name != "numpy":
        raise
    np = None


@unittest.skipIf(np is None, "NumPy is required for RMSD calculation tests")
class RMSDTests(unittest.TestCase):
    def test_rigid_motion_is_removed_but_reflection_is_not(self):
        reference = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
        rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        mobile = reference @ rotation + np.array([4., -2., 5.])
        self.assertLess(rmsd.kabsch_rmsd(reference, mobile), 1e-12)
        mirrored = reference * np.array([-1., 1., 1.])
        self.assertGreater(rmsd.kabsch_rmsd(reference, mirrored), 0.1)

    def test_hydrogen_motion_is_excluded_and_atom_mapping_is_checked(self):
        reference = rmsd.XYZGeometry(("C", "O", "H"), np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]]))
        mobile = rmsd.XYZGeometry(reference.symbols, reference.coordinates.copy())
        mobile.coordinates[2] = [100., 200., 300.]
        ref_heavy, mobile_heavy, count = rmsd.heavy_atom_coordinates(reference, mobile)
        self.assertEqual(count, 2)
        self.assertLess(rmsd.kabsch_rmsd(ref_heavy, mobile_heavy), 1e-12)
        for symbols in (("O", "C", "H"), ("C", "O")):
            with self.subTest(symbols=symbols), self.assertRaises(ValueError):
                rmsd.validate_atom_mapping(reference, rmsd.XYZGeometry(symbols, mobile.coordinates))

    def test_multiframe_and_nonfinite_inputs_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.xyz"
            for content in ("1\nfirst\nC 0 0 0\n1\nsecond\nC 0 0 0\n", "1\ncomment\nC nan 0 0\n"):
                path.write_text(content)
                with self.subTest(content=content), self.assertRaises(ValueError):
                    rmsd.read_xyz(path)

    def test_only_the_current_exact_filename_is_selected(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            (folder / "ac16_major_optimized_structure.xyz").write_text("legacy")
            (folder / "unrelated.xyz").write_text("unrelated")
            self.assertIsNone(rmsd.choose_xyz(folder, "ac16_major"))
            canonical = folder / "ac16_major.xyz"
            canonical.write_text("current")
            self.assertEqual(rmsd.choose_xyz(folder, "ac16_major"), canonical)

    def test_cli_handles_current_mace_names_and_reports_missing_structures(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            content = "3\ncomment\nC 0 0 0\nO 1 0 0\nH 0 1 0\n"
            for method in (rmsd.DEFAULT_REFERENCE, "MACE-OMOL"):
                for structure in ("ac16_major", "ac16_minor"):
                    folder = root / method / structure
                    folder.mkdir(parents=True)
                    (folder / f"{structure}.xyz").write_text(content)
            output = root / "rmsd.csv"
            command = [sys.executable, str(Path(rmsd.__file__).resolve()), "--root", str(root),
                       "--output", str(output), "--strict"]
            complete = subprocess.run(command, cwd=directory, capture_output=True, text=True)
            self.assertEqual(complete.returncode, 0, complete.stderr)
            with output.open(encoding="utf-8-sig", newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 2)
            self.assertTrue(all(r["status"] == "ok" and r["heavy_atom_count"] == "2" for r in rows))
            self.assertTrue(all(float(r["rmsd_angstrom"]) == 0 for r in rows))
            (root / "MACE-OMOL/ac16_minor/ac16_minor.xyz").unlink()
            missing = subprocess.run(command, cwd=directory, capture_output=True, text=True)
            self.assertEqual(missing.returncode, 1, missing.stderr)
            with output.open(encoding="utf-8-sig", newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([r["structure_name"] for r in rows if r["status"] == "error"], ["ac16_minor"])


if __name__ == "__main__":
    unittest.main()
