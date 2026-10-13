import copy
import tempfile
import unittest
from pathlib import Path

from scripts.validate_data import (REFERENCE, ENERGY_METHODS, check_bond, check_bond_coverage,
                                   check_energy_orientation, check_links, number, read_table, read_xyz)


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

    def test_bond_coverage_allows_only_declared_exclusions(self):
        structures = {f"s{i}" for i in range(98)} | {"ja1g_major", "ja1g_minor"}
        bonds = {(f"s{i}", 1, 2) for i in range(98)} | {(f"s{i}", 2, 3) for i in range(68)}
        failures = {("s0", "other")}
        available = {b for b in bonds if b[0] != "s0"}
        check_bond_coverage({REFERENCE: bonds, "other": available}, structures, failures)
        for wrong in (bonds, available - {("s1", 1, 2)}):
            with self.subTest(wrong=len(wrong)), self.assertRaisesRegex(ValueError, "coverage differs"):
                check_bond_coverage({REFERENCE: bonds, "other": wrong}, structures, failures)
        with self.assertRaisesRegex(ValueError, "98 structures"):
            check_bond_coverage({REFERENCE: bonds}, structures | {"extra"})
        with self.assertRaisesRegex(ValueError, "98 structures"):
            check_bond_coverage({REFERENCE: bonds | {("ja1g_minor", 1, 3)}}, structures)

    def test_nullable_columns_still_reject_missing_fields(self):
        self.assertEqual(read_table(self.write("valid.csv", "id,value\na,\n"),
                                    ["id", "value"], nullable=("value",)), [{"id": "a", "value": ""}])
        for text in ("id,value\na\n", "id,value\n,1\n", "id,value\na,1,extra\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                read_table(self.write("invalid.csv", text), ["id", "value"], nullable=("value",))

    def test_required_schema_allows_annotations_but_rejects_missing_or_duplicate_columns(self):
        rows = read_table(self.write("annotated.csv", "id,value,source\na,1,\n"),
                          ["id", "value"], allow_extra=True)
        self.assertEqual(rows[0]["source"], "")
        for text in ("id,source\na,1\n", "id,value,value\na,1,2\n", "id,value\na,\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                read_table(self.write("bad.csv", text), ["id", "value"], allow_extra=True)

    def orientations(self):
        energies, rows = {}, []
        for q in ("delta_delta_e", "delta_delta_g"):
            energies[q] = {}
            for case in ("example", "an82"):
                factor = -1 if case == "example" else 1
                gap = 2 if q == "delta_delta_e" else -0.5
                current = {f"dft_reference_{q}_kcal_mol": str(gap)}
                for model, prefix in ENERGY_METHODS.items():
                    value = "" if case == "an82" and model == "gxtb" else str(gap + 0.1)
                    current[f"{prefix}_{q}_kcal_mol"] = value
                    rows.append(dict(quantity=q, benchmark_id=case, model_id=model,
                                     source_reference_gap_kcal_mol=str(gap * factor),
                                     source_predicted_gap_kcal_mol=str(float(value) * factor) if value else "",
                                     orientation_factor=str(factor)))
                energies[q][case] = current
        return rows, energies, {"example", "an82"}

    def test_fixed_roles_allow_negative_g_and_declared_failed_values(self):
        check_energy_orientation(*self.orientations())

    def test_orientation_rejects_g_role_flip_and_missing_or_duplicate_records(self):
        rows, energies, cases = self.orientations()
        bad = copy.deepcopy(rows)
        for row in bad:
            if row["quantity"] == "delta_delta_g" and row["benchmark_id"] == "example":
                row["orientation_factor"] = "1"
        for wrong in (bad, rows[:-1], rows + [rows[0]]):
            with self.subTest(count=len(wrong)), self.assertRaises(ValueError):
                check_energy_orientation(wrong, energies, cases)

    def test_orientation_rejects_undeclared_missing_values_and_imputed_failures(self):
        rows, energies, cases = self.orientations()
        for case, model, value in (("example", "esen", ""), ("an82", "gxtb", "0")):
            bad = copy.deepcopy(energies)
            bad["delta_delta_e"][case][f"{ENERGY_METHODS[model]}_delta_delta_e_kcal_mol"] = value
            with self.subTest(case=case), self.assertRaises(ValueError):
                check_energy_orientation(rows, bad, cases)

    def test_local_links_allow_fragments_and_parentheses_but_reject_missing_files(self):
        self.write("method(BJ).xyz", "content")
        self.write("README.md", "[data](method(BJ).xyz) [here](#section) [web](https://example.com)\n")
        check_links(self.root)
        self.write("README.md", "[missing](absent.csv)\n")
        with self.assertRaisesRegex(ValueError, "broken link"):
            check_links(self.root)


if __name__ == "__main__":
    unittest.main()
