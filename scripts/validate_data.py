"""Validate the supplied SELECT50 snapshot without external dependencies."""

import csv
import math
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
METHODS = {
    "esen": "eSEN",
    "gxtb": "g-xTB",
    "mace_omol": "MACE-OMOL",
    "pbe0_d3bj_svp": "PBE0-D3(BJ)-def2-SVP",
    "r2scan_3c": "r2SCAN-3c",
    "uma_m": "UMA-M-1.1",
    "uma_s": "UMA-S-1.2",
    "wb97x_d3bj_svp": "wB97X-D3(BJ)-def2-SVP",
    "dft_geometry_reference": "wB97X-D3(BJ)-def2-TZVP",
}
REFERENCE = "dft_geometry_reference"
# Author confirmed 2026-10-12: preserve the current SI/table exclusion of both.
BONDLESS_STRUCTURES = {"ja1g_major", "ja1g_minor"}
GEOMETRY_FAILURES = {(s, m) for s in ("an82_major", "an82_minor")
                     for m in ("gxtb", "r2scan_3c")}
ENERGY_FAILURES = {(q, "an82", "gxtb") for q in ("delta_delta_e", "delta_delta_g")}
ENERGY_METHODS = {"esen": "esen_medium", "mace_omol": "mace_omol",
                  "uma_m": "uma_m_1_1", "uma_s": "uma_s_1_2",
                  "gxtb": "g_xtb", "specorr": "specorr"}
ELEMENTS = set("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def number(value):
    result = float(value)
    require(math.isfinite(result), f"Non-finite number: {value}")
    return result


def read_xyz(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    frames = []
    offset = 0
    while offset < len(lines):
        if not lines[offset].strip():
            offset += 1
            continue
        count = int(lines[offset])
        require(count > 0 and offset + count + 2 <= len(lines), f"{path}: incomplete XYZ frame")
        atoms = []
        for line in lines[offset + 2:offset + count + 2]:
            fields = line.split()
            require(len(fields) >= 4 and fields[0] in ELEMENTS, f"{path}: invalid atom: {line}")
            atoms.append((fields[0], tuple(number(v) for v in fields[1:4])))
        frames.append(atoms)
        offset += count + 2
    require(frames, f"{path}: empty XYZ file")
    return frames


def read_table(path, columns, nullable=(), allow_extra=False):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames is not None and len(reader.fieldnames) == len(set(reader.fieldnames))
                and (set(columns) <= set(reader.fieldnames) if allow_extra else reader.fieldnames == columns),
                f"{path}: unexpected CSV columns")
        rows = list(reader)
    require(rows, f"{path}: empty table")
    for row in rows:
        require(None not in row and all(v is not None and (v.strip() or k in nullable or (allow_extra and k not in columns)) for k, v in row.items()),
                f"{path}: missing or extra CSV fields")
    return rows


def check_bond(row, atoms):
    i = int(row["atom_i_1based"]) - 1
    j = int(row["atom_j_1based"]) - 1
    require(0 <= i < len(atoms) and 0 <= j < len(atoms) and i != j,
            f"Invalid bond indices: {row}")
    require((atoms[i][0], atoms[j][0]) == (row["atom_i_symbol"], row["atom_j_symbol"]),
            f"Bond atom symbols disagree with XYZ: {row}")
    distance = number(row["distance_angstrom"])
    require(distance > 0, f"Invalid bond distance: {row}")
    actual = math.dist(atoms[i][1], atoms[j][1])
    require(abs(actual - distance) <= 0.00002,
            f"Bond distance disagrees with XYZ ({actual:.6f}): {row}")


def check_links(root):
    for path in root.rglob("*.md"):
        if any(part.startswith(".") for part in path.relative_to(root).parts):
            continue
        # Current docs use inline links; allow one level of parentheses in paths.
        for target in re.findall(r"\[[^\]]*\]\(([^()\s]*(?:\([^()]*\)[^()\s]*)*)\)", path.read_text()):
            url = urlsplit(target)
            if url.scheme or url.netloc or not url.path:
                continue
            require((path.parent / unquote(url.path)).exists(), f"{path}: broken link: {target}")


def check_bond_coverage(bond_sets, structures, failures=()):
    reference_bonds = bond_sets[REFERENCE]
    require(len(reference_bonds) == 166 and
            {s for s, _, _ in reference_bonds} == structures - BONDLESS_STRUCTURES,
            "Expected 166 reference bonds covering 98 structures; both ja1g structures are excluded")
    for method, bonds in bond_sets.items():
        expected = {b for b in reference_bonds if (b[0], method) not in failures}
        require(bonds == expected, f"{method}: bond coverage differs from reference minus declared failures")


def check_energy_orientation(rows, energies, cases):
    expected = {(q, c, m) for q in energies for c in cases for m in ENERGY_METHODS}
    seen, factors = set(), {}
    for row in rows:
        q, c, m = row["quantity"], row["benchmark_id"], row["model_id"]
        key = q, c, m
        require(key in expected and key not in seen, f"Unexpected/duplicate orientation: {key}")
        seen.add(key)
        factor = number(row["orientation_factor"])
        require(factor in (-1, 1), f"Invalid orientation: {key}")
        require(c not in factors or factors[c] == factor, f"E/G roles differ: {c}")
        factors[c] = factor
        current = energies[q][c]
        for source, column in (("source_reference_gap_kcal_mol", f"dft_reference_{q}_kcal_mol"),
                               ("source_predicted_gap_kcal_mol", f"{ENERGY_METHODS[m]}_{q}_kcal_mol")):
            if key in ENERGY_FAILURES and source == "source_predicted_gap_kcal_mol":
                require(row[source] == current[column] == "", f"Failed energy must remain blank: {key}")
            else:
                require(abs(number(row[source]) * factor - number(current[column])) <= 1e-6,
                        f"Source orientation disagrees with canonical energy: {key}/{column}")
    require(seen == expected, "Missing energy orientation records")


def check_energy_summary(rows, energies):
    expected = {(q, m) for q in energies for m in ENERGY_METHODS}
    seen = set()
    for row in rows:
        q, m = row["quantity"], row["model_id"]
        require((q, m) in expected and (q, m) not in seen, "Unexpected/duplicate energy summary")
        seen.add((q, m))
        pairs = [(number(r[f"{ENERGY_METHODS[m]}_{q}_kcal_mol"]),
                  number(r[f"dft_reference_{q}_kcal_mol"])) for r in energies[q].values()
                 if r[f"{ENERGY_METHODS[m]}_{q}_kcal_mol"]]
        errors = [p - r for p, r in pairs]
        n = len(errors)
        correct = sum((p > 0) == (r > 0) for p, r in pairs)
        values = {"n_pairs_in_category": len(energies[q]), "n_pairs_evaluated": n,
                  "mae_kcal_mol": sum(abs(e) for e in errors) / n,
                  "rmse_kcal_mol": math.sqrt(sum(e * e for e in errors) / n),
                  "mean_signed_error_kcal_mol": sum(errors) / n,
                  "n_correct": correct, "accuracy_percent": 100 * correct / n}
        for field, value in values.items():
            require(abs(number(row[field]) - value) <= 1e-5, f"Energy summary differs: {q}/{m}/{field}")
    require(seen == expected, "Missing energy summaries")


def validate(root=ROOT):
    reference_root = root / "results/reference"
    failure_rows = read_table(reference_root / "failed_calculations.csv",
                             ["quantity", "benchmark_id", "structure_id", "model_id", "status", "reason", "evidence"],
                             nullable=("structure_id",))
    failures = {(r["structure_id"], r["model_id"]) for r in failure_rows if r["quantity"] == "geometry"}
    failed_energies = ENERGY_FAILURES
    require(len(failure_rows) == 4 and failures == GEOMETRY_FAILURES
            and all(r["status"] == "calculation_failed" and r["benchmark_id"] == "an82" for r in failure_rows),
            "Failure ledger differs from the four documented an82 geometry exclusions")
    # Required analysis columns are fixed here; additional source annotations are allowed.
    energy_metadata = ["benchmark_id", "description", "molecular_formula", "selectivity_class", "doi",
                       "charge_e", "spin_multiplicity", "total_atom_count", "heavy_atom_count"]
    bond_columns = ["benchmark_id", "structure_id", "model_id", "reference_model_id", "atom_i_1based",
                    "atom_j_1based", "atom_lo_1based", "atom_hi_1based", "bond_id", "atom_i_symbol",
                    "atom_j_symbol", "distance_angstrom"]
    required_columns = {
        "geometry/heavy_atom_rmsd.csv": ["structure_id", "model_id", "reference_model_id", "benchmark_id",
                                        "heavy_atom_count", "rmsd_angstrom"],
        "geometry/ts_bond_distances.csv": bond_columns,
        "geometry/ts_bond_errors.csv": bond_columns + ["reference_distance_angstrom", "signed_error_angstrom", "abs_error_angstrom"],
        "geometry/heavy_atom_rmsd_method_summary.csv": ["model_id", "n_structures", "mean_rmsd_angstrom", "median_rmsd_angstrom"],
        "geometry/ts_bond_method_summary.csv": ["model_id", "n_bonds", "n_structures", "bond_error_rmse_angstrom", "mae_angstrom"],
        "geometry/ts_bond_category_summary.csv": ["table_id", "category_field", "category", "model_id", "n_systems", "n_bonds", "mae_angstrom"],
        "energetics/energy_method_summary.csv": ["quantity", "model_id", "n_pairs_in_category", "n_pairs_evaluated",
                                                "mae_kcal_mol", "rmse_kcal_mol", "mean_signed_error_kcal_mol", "n_correct", "accuracy_percent"],
        "energetics/energy_category_mae.csv": ["quantity", "table_id", "category_field", "category", "model_id",
                                              "n_pairs_in_category", "n_pairs_evaluated", "mae_kcal_mol"],
        "energetics/energy_ordering_accuracy.csv": ["quantity", "table_id", "model_id", "gap_bin", "n_pairs_in_bin",
                                                   "n_pairs_evaluated", "n_correct", "accuracy_percent"],
    }
    nullable = ["doi"]
    for q in ("delta_delta_e", "delta_delta_g"):
        required_columns[f"energetics/{q}.csv"] = energy_metadata + [f"{m}_{q}_kcal_mol" for m in ["dft_reference", *ENERGY_METHODS.values()]]
        nullable.append(f"g_xtb_{q}_kcal_mol")
    tables = {f"results/{path}": read_table(root / "results" / path, columns, nullable, allow_extra=True)
              for path, columns in required_columns.items()}
    xyz = {p: read_xyz(p) for p in root.rglob("*.xyz")}
    inputs = reference_root / "starting_structures"
    refs = {p.stem: frames[0] for p, frames in xyz.items() if p.parent == inputs}
    require(len(refs) == 100, "Expected 100 reference structures")
    require(all(re.fullmatch(r"[a-z0-9]+_(major|minor)", s) for s in refs), "Invalid structure ID")
    cases = {s.rsplit("_", 1)[0] for s in refs}
    require(len(cases) == 50 and "an82" in cases and "ja30" not in cases and
            set(refs) == {f"{c}_{role}" for c in cases for role in ("major", "minor")},
            "Expected 50 complete current major/minor pairs including an82, excluding ja30")
    metadata = read_table(reference_root / "case_metadata.csv",
                          ["benchmark_id", "description", "molecular_formula", "selectivity_class", "charge_e",
                           "spin_multiplicity", "doi", "metal_elements", "metal_group", "metal_containing"],
                          nullable=("doi", "metal_elements", "metal_group"))
    require(len(metadata) == 50 and {r["benchmark_id"] for r in metadata} == cases, "Case metadata coverage differs")
    metadata = {r["benchmark_id"]: r for r in metadata}
    geometries = {}
    folders = {p.name for p in (root / "optimized-geometries").iterdir() if p.is_dir()}
    require(set(METHODS.values()) == folders, "Optimized method folders differ from registered methods")
    for model, folder in METHODS.items():
        method_root = root / "optimized-geometries" / folder
        method = {p.stem: frames[0] for p, frames in xyz.items() if method_root in p.parents}
        expected = {s for s in refs if (s, model) not in failures}
        require(set(method) == expected, f"{folder}: structure IDs differ from reference minus declared failures")
        paths = {method_root / s / f"{s}.xyz" for s in expected}
        require({p for p in xyz if method_root in p.parents} == paths, f"{folder}: unexpected or duplicate XYZ paths")
        for s, atoms in method.items():
            require([a[0] for a in atoms] == [a[0] for a in refs[s]], f"{folder}/{s}: atom ordering differs")
        geometries[model] = method
    for p, frames in xyz.items():
        if inputs in p.parents or (root / "optimized-geometries") in p.parents:
            require(len(frames) == 1, f"{p}: expected one XYZ frame")
    for s in refs:
        copy = root / "optimized-geometries" / METHODS[REFERENCE] / s / f"{s}.xyz"
        require((inputs / f"{s}.xyz").read_bytes() == copy.read_bytes(), f"{s}: reference copies differ")
    for relative_path, count in (("cycloaddition/cycloaddition_structures.xyz", 10),
                                 ("palladium_cyclization/structures/uma/all_21_frames.xyz", 21)):
        p = root / "reaction_landscapes" / relative_path
        require(p in xyz and len(xyz[p]) == count, f"{relative_path}: unexpected landscape frame count")
    energies = {}
    for q in ("delta_delta_e", "delta_delta_g"):
        rows = tables[f"results/energetics/{q}.csv"]
        require(len(rows) == 50 and {r["benchmark_id"] for r in rows} == cases, f"{q}: missing or duplicate cases")
        energies[q] = {r["benchmark_id"]: r for r in rows}
        for row in rows:
            c = row["benchmark_id"]
            ref = number(row[f"dft_reference_{q}_kcal_mol"])
            require(q != "delta_delta_e" or ref > 0, f"{c}: reference E roles reversed")
            for m, prefix in ENERGY_METHODS.items():
                value = row[f"{prefix}_{q}_kcal_mol"]
                if (q, c, m) in failed_energies:
                    require(value == "", f"{q}/{c}/{m}: declared failed energy must be blank")
                else:
                    number(value)
            for field in ("description", "molecular_formula", "selectivity_class", "charge_e", "spin_multiplicity", "doi"):
                require(row[field] == metadata[c][field], f"{q}/{c}: metadata differs for {field}")
            for role in ("major", "minor"):
                atoms = refs[f"{c}_{role}"]
                require(int(row["total_atom_count"]) == len(atoms) and
                        int(row["heavy_atom_count"]) == sum(a[0] != "H" for a in atoms), f"{c}: metadata atom counts differ")
    orientation = read_table(reference_root / "energy_orientation.csv",
                              ["benchmark_id", "quantity", "model_id", "source_reference_gap_kcal_mol",
                               "source_predicted_gap_kcal_mol", "orientation_factor"],
                              nullable=("source_predicted_gap_kcal_mol",))
    check_energy_orientation(orientation, energies, cases)
    check_energy_summary(tables["results/energetics/energy_method_summary.csv"], energies)
    rmsd = tables["results/geometry/heavy_atom_rmsd.csv"]
    expected = {(s, t) for s in refs for t in METHODS if t != REFERENCE} - failures
    require(len(rmsd) == len(expected) and {(r["structure_id"], r["model_id"]) for r in rmsd} == expected,
            "RMSD table has missing, duplicate, or unknown structure/method pairs")
    for row in rmsd:
        atoms = refs[row["structure_id"]]
        require(row["reference_model_id"] == REFERENCE and row["benchmark_id"] == row["structure_id"].rsplit("_", 1)[0],
                f"RMSD identity mismatch: {row}")
        require(int(row["heavy_atom_count"]) == sum(a[0] != "H" for a in atoms), f"Wrong heavy-atom count: {row}")
        require(number(row["rmsd_angstrom"]) >= 0, f"Negative RMSD: {row}")
    bonds = tables["results/geometry/ts_bond_distances.csv"]
    measured, bond_sets = {}, {t: set() for t in METHODS}
    for row in bonds:
        s, m = row["structure_id"], row["model_id"]
        i, j = int(row["atom_i_1based"]), int(row["atom_j_1based"])
        key = s, m, min(i, j), max(i, j)
        require(s in refs and m in METHODS and s in geometries[m], f"Unknown bond structure/method: {key}")
        require(key not in measured, f"Duplicate bond: {key}")
        require(row["reference_model_id"] == REFERENCE and row["benchmark_id"] == s.rsplit("_", 1)[0]
                and (int(row["atom_lo_1based"]), int(row["atom_hi_1based"])) == key[2:], f"Bond identity mismatch: {key}")
        check_bond(row, geometries[m][s])
        measured[key] = number(row["distance_angstrom"])
        bond_sets[m].add((s, key[2], key[3]))
    check_bond_coverage(bond_sets, set(refs), failures)
    seen = set()
    for row in tables["results/geometry/ts_bond_errors.csv"]:
        s, m = row["structure_id"], row["model_id"]
        i, j = sorted((int(row["atom_i_1based"]), int(row["atom_j_1based"])))
        key = s, m, i, j
        require(m != REFERENCE and key in measured and key not in seen, f"Unexpected/duplicate bond error: {key}")
        seen.add(key)
        reference = measured[(s, REFERENCE, i, j)]
        error = measured[key] - reference
        for field, actual in (("distance_angstrom", measured[key]), ("reference_distance_angstrom", reference),
                              ("signed_error_angstrom", error), ("abs_error_angstrom", abs(error))):
            require(abs(number(row[field]) - actual) <= 2e-6, f"Bond error differs: {key}/{field}")
    require(seen == {k for k in measured if k[1] != REFERENCE}, "Missing bond error rows")
    check_links(root)
    print(f"PASS: {len(cases)} cases, {len(xyz)} XYZ files, {sum(map(len, xyz.values()))} frames, "
          f"{len(folders)} method sets, {len(rmsd)} RMSD rows, {len(bonds)} verified bond rows; local links valid.")
    print("PASS: required schemas, metadata counts, 600 orientations, energy summaries and 1322 bond errors.")
    print("NOTE: two failed energy values and four failed geometries are declared for an82; both ja1g structures are excluded from bond analysis.")


if __name__ == "__main__":
    try:
        validate()
    except (ValueError, OSError, KeyError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        sys.exit(1)
