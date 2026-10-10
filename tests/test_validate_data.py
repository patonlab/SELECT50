import tempfile
import unittest
from pathlib import Path

from scripts.validate_data import REFERENCE, check_bond, check_bond_coverage, check_links, number, read_table, read_xyz


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write(self, name, text):
        path = self.root / name
        path.write_text(text)
        return path

    def test_xyz_reads_multiple_frames_with_blank_comments(self):
        path = self.write("frames.xyz", "2\n\nC 0 0 0\nH 0 0 1\n\n1\nsecond\nHe 1 2 3\n")
        frames = read_xyz(path)
        self.assertEqual([len(f) for f in frames], [2, 1])
        self.assertEqual(frames[1][0], ("He", (1.0, 2.0, 3.0)))

    def test_xyz_rejects_truncation_unknown_elements_and_nonfinite_coordinates(self):
        for text in ("2\ncomment\nH 0 0 0\n", "1\ncomment\nXx 0 0 0\n", "1\ncomment\nH nan 0 0\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                read_xyz(self.write("bad.xyz", text))

    def test_csv_rejects_ragged_rows_and_changed_schema(self):
        for text in ("id,value\na\n", "id,value\na,1,extra\n", "id,value\na,\n", "value,id\n1,a\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                read_table(self.write("bad.csv", text), ["id", "value"])
        self.assertEqual(read_table(self.write("good.csv", "id,value\na,1\n"), ["id", "value"]),
                         [{"id": "a", "value": "1"}])

    def test_bonds_use_one_based_indices_and_validate_coordinates(self):
        atoms = [("C", (0, 0, 0)), ("H", (0, 0, 1))]
        row = {"atom_i_1based": "1", "atom_j_1based": "2", "atom_i_symbol": "C",
               "atom_j_symbol": "H", "distance_angstrom": "1.000001"}
        check_bond(row, atoms)
        for changes in ({"atom_i_1based": "0"}, {"atom_j_1based": "3"},
                        {"atom_i_symbol": "N"}, {"distance_angstrom": "1.1"},
                        {"distance_angstrom": "nan"}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                check_bond(row | changes, atoms)

    def test_nonfinite_numbers_are_rejected(self):
        for value in ("nan", "inf", "-inf"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                number(value)

    def test_bond_coverage_rejects_missing_structures_and_method_rows(self):
        structures = {f"s{i}" for i in range(98)} | {"ja1g_major", "ja1g_minor"}
        bonds = {(f"s{i}", 1, 2) for i in range(98)} | {(f"s{i}", 2, 3) for i in range(69)}
        check_bond_coverage({REFERENCE: bonds, "other": bonds}, structures)
        with self.assertRaisesRegex(ValueError, "differs between methods"):
            check_bond_coverage({REFERENCE: bonds, "other": bonds - {("s0", 1, 2)}}, structures)
        with self.assertRaisesRegex(ValueError, "98 structures"):
            check_bond_coverage({REFERENCE: bonds}, structures | {"extra"})

    def test_local_links_allow_fragments_and_parentheses_but_reject_missing_files(self):
        self.write("method(BJ).xyz", "content")
        self.write("README.md", "[data](method(BJ).xyz) [here](#section) [web](https://example.com)\n")
        check_links(self.root)
        self.write("README.md", "[missing](absent.csv)\n")
        with self.assertRaisesRegex(ValueError, "broken link"):
            check_links(self.root)


if __name__ == "__main__":
    unittest.main()
