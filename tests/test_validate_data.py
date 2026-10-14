import copy
import csv
import shutil
import tempfile
import unittest
from pathlib import Path

from scripts.validate_data import (REFERENCE, ENERGY_METHODS, check_bond, check_bond_coverage,
                                   check_energy_categories, check_energy_ordering, check_raw_energy_table,
                                   check_raw_energies, check_energy_orientation, check_links, number, read_table, read_xyz)


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

    def test_raw_energies_reject_duplicates_missing_records_and_sign_changes(self):
        rows = [dict(structure_id="example_major", E_hartree="-1.000"),
                dict(structure_id="example_minor", E_hartree="-1.001")]
        energies = {"delta_delta_e": {"example": {"esen_medium_delta_delta_e_kcal_mol": "-0.6275094740631"}}}
        def check(raw):
            return check_raw_energy_table(raw, energies, "delta_delta_e", "esen",
                                          "structure_id", "E_hartree", "example.csv")
        self.assertEqual(check(rows)[1], 1)
        for bad in (rows + [rows[0]], rows[:1],
                    [rows[0], rows[1] | {"E_hartree": "-0.999"}],
                    [rows[0], rows[1] | {"E_hartree": "nan"}],
                    [rows[0], rows[1] | {"structure_id": "example_min"}]):
            with self.subTest(rows=bad), self.assertRaises(ValueError):
                check(bad)
        energies["delta_delta_e"]["example"]["esen_medium_delta_delta_e_kcal_mol"] = "-0.627"
        check(rows)

    def test_raw_exclusions_allow_an82_without_hiding_other_missing_pairs(self):
        energies = {"delta_delta_g": {
            "example": {"g_xtb_delta_delta_g_kcal_mol": "0.6275094740631"},
            "an82": {"g_xtb_delta_delta_g_kcal_mol": ""}}}
        rows = [dict(structure_id="example_major", G_hartree="-1.001"),
                dict(structure_id="example_minor", G_hartree="-1.000"),
                dict(structure_id="an82_minor", G_hartree="-2")]
        args = (energies, "delta_delta_g", "gxtb", "structure_id", "G_hartree", "gxtb.csv")
        self.assertEqual(check_raw_energy_table(rows, *args)[1], 1)
        with self.assertRaisesRegex(ValueError, "missing raw"):
            check_raw_energy_table(rows[1:], *args)

    def energy_tables(self):
        base = Path(__file__).resolve().parents[1] / "results/energetics"
        def load(name):
            with (base / name).open(encoding="utf-8-sig", newline="") as handle:
                return list(csv.DictReader(handle))
        energies = {q: {r['benchmark_id']: r for r in load(q + '.csv')}
                    for q in ('delta_delta_e', 'delta_delta_g')}
        return load, energies

    def test_category_summaries_reject_wrong_values_and_incomplete_coverage(self):
        load, energies = self.energy_tables()
        rows = load("energy_category_mae.csv")
        check_energy_categories(rows, energies)
        wrong = copy.deepcopy(rows)
        wrong[0]['mae_kcal_mol'] = str(float(wrong[0]['mae_kcal_mol']) + .01)
        for bad in (wrong, rows[:-1], rows + [rows[0]]):
            with self.subTest(count=len(bad)), self.assertRaises(ValueError):
                check_energy_categories(bad, energies)

    def test_ordering_summaries_reject_wrong_counts_and_incomplete_coverage(self):
        load, energies = self.energy_tables()
        rows = load("energy_ordering_accuracy.csv")
        check_energy_ordering(rows, energies)
        wrong = copy.deepcopy(rows)
        wrong[0]['n_correct'] = str(int(wrong[0]['n_correct']) - 1)
        for bad in (wrong, rows[:-1], rows + [rows[0]]):
            with self.subTest(count=len(bad)), self.assertRaises(ValueError):
                check_energy_ordering(bad, energies)

    def test_ordering_bin_boundaries_use_absolute_reference_gaps(self):
        q = 'delta_delta_e'
        energies = {q: {}}
        for i, (reference, prediction) in enumerate([(-.49, -.4), (.5, -.3), (-1, 1), (1.01, 2)]):
            energies[q][str(i)] = {f'dft_reference_{q}_kcal_mol': str(reference)} | {
                f'{prefix}_{q}_kcal_mol': str(prediction) for prefix in ENERGY_METHODS.values()}
        rows = [dict(quantity=q, table_id='S13', model_id=model, gap_bin=bin_id,
                     n_pairs_in_bin=str(n), n_pairs_evaluated=str(n), n_correct=str(correct),
                     accuracy_percent=str(100 * correct / n))
                for model in ENERGY_METHODS if model != 'specorr'
                for bin_id, n, correct in [('all', 4, 2), ('lt_0p5', 1, 1),
                                          ('0p5_to_1p0_inclusive', 2, 0), ('gt_1p0', 1, 1)]]
        check_energy_ordering(rows, energies)

    def test_raw_file_validation_checks_mace_thermochemistry_and_electronic_consistency(self):
        _, energies = self.energy_tables()
        source = Path(__file__).resolve().parents[1] / 'results/energetics'
        target = self.root / 'results/energetics'
        target.mkdir(parents=True)
        for path in source.glob('*_raw_*.csv'):
            shutil.copy2(path, target)
        self.assertEqual(check_raw_energies(self.root, energies), 598)
        path = target / 'mace_omol_raw_thermochemistry.csv'
        with path.open(newline='') as handle:
            rows = list(csv.DictReader(handle))
        bad_correction = copy.deepcopy(rows)
        bad_correction[0]['thermal_correction_hartree'] = 'nan'
        broken_sum = copy.deepcopy(rows)
        broken_sum[0]['thermal_correction_hartree'] = str(float(rows[0]['thermal_correction_hartree']) + .01)
        common_offset = copy.deepcopy(rows)
        for row in common_offset:
            for field in ('electronic_energy_hartree', 'free_energy_hartree'):
                row[field] = str(float(row[field]) + .01)
        for bad in (bad_correction, broken_sum, common_offset):
            with self.subTest(first=bad[0]), self.assertRaises(ValueError):
                with path.open('w', newline='') as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(bad[0]))
                    writer.writeheader()
                    writer.writerows(bad)
                check_raw_energies(self.root, energies)


if __name__ == "__main__":
    unittest.main()
