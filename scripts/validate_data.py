"""Validate the supplied SELECT50 snapshot without external dependencies."""

import csv
import math
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
METHODS = {
    "eSEN-medium": "eSEN",
    "g-xTB": "g-xTB",
    "MACE-OMOL": "MACE-OMOL",
    "PBE0-D3(BJ)/def2-SVP": "PBE0-D3(BJ)-def2-SVP",
    "r2SCAN-3c": "r2SCAN-3c",
    "UMA-M-1.2": "UMA-M",
    "UMA-S-1.2": "UMA-S",
    "wB97X-D3(BJ)/def2-SVP": "wB97X-D3(BJ)-def2-SVP",
    "wB97X-D3(BJ)/def2-TZVP": "wB97X-D3(BJ)-def2-TZVP",
}
REFERENCE = "wB97X-D3(BJ)/def2-TZVP"
# The supplied render identifies ja1g_major as rotation, with no assigned TS bond.
BONDLESS_STRUCTURES = {"ja1g_major"}
ELEMENTS = set("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split())
ENERGY_METHODS = ("dft", "esen", "mace_omol", "uma_m", "uma_s", "gxtb", "spe_corr")


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


def read_table(path, columns):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames == columns, f"{path}: unexpected CSV columns")
        rows = list(reader)
    require(rows, f"{path}: empty table")
    for row in rows:
        require(None not in row and all(v is not None and v.strip() for v in row.values()),
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


def check_bond_coverage(bond_sets, structures):
    reference_bonds = bond_sets[REFERENCE]
    require(len(reference_bonds) == 168 and {s for s, _, _ in reference_bonds} == structures - BONDLESS_STRUCTURES,
            "Expected 168 reference bonds covering 99 structures; ja1g_major is a rotational TS without an assigned bond")
    require(("ja1g_minor", 1, 3) in reference_bonds, "ja1g_minor requires the assigned P1–C3 breaking bond")
    require(all(b == reference_bonds for b in bond_sets.values()), "Bond coverage differs between methods")


def validate(root=ROOT):
    xyz = {p: read_xyz(p) for p in root.rglob("*.xyz")}
    inputs = root / "results/reference/starting_structures"
    refs = {p.stem: frames[0] for p, frames in xyz.items() if p.parent == inputs}
    require(len(refs) == 100, "Expected 100 reference structures")
    require(all(re.fullmatch(r"[a-z0-9]+_(major|minor)", s) for s in refs), "Invalid structure ID")
    cases = {s.rsplit("_", 1)[0] for s in refs}
    require(len(cases) == 50 and set(refs) == {f"{c}_{role}" for c in cases for role in ("major", "minor")},
            "Expected 50 complete major/minor pairs")
    geometries = {}
    folders = {p.name for p in (root / "optimized-geometries").iterdir() if p.is_dir()}
    require(set(METHODS.values()) <= folders, "Missing optimized method folders")
    for folder in sorted(folders):
        method_root = root / "optimized-geometries" / folder
        method = {p.stem: frames[0] for p, frames in xyz.items() if method_root in p.parents}
        require(set(method) == set(refs), f"{folder}: structure IDs differ from reference")
        paths = {method_root / s / f"{s}.xyz" for s in refs}
        actual_paths = {p for p in xyz if method_root in p.parents}
        require(actual_paths == paths, f"{folder}: unexpected or duplicate XYZ paths")
        for s, atoms in method.items():
            require([a[0] for a in atoms] == [a[0] for a in refs[s]], f"{folder}/{s}: atom ordering differs")
        geometries[folder] = method
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

    for kind in ("e", "g"):
        path = root / f"results/energetics/delta_{kind}.csv"
        columns = ["benchmark_id"] + [f"{m}_delta_{kind}_minor_minus_major (kcal/mol)" for m in ENERGY_METHODS]
        rows = read_table(path, columns)
        require(len(rows) == 50 and {r["benchmark_id"] for r in rows} == cases,
                f"{path}: missing or duplicate cases")
        for row in rows:
            for col in columns[1:]:
                number(row[col])
            if kind == "e":
                require(number(row[columns[1]]) > 0, f"{path}: reference ordering is not minor-minus-major")

    rmsd = read_table(root / "results/geometry/heavy_atom_rmsd.csv",
                      ["structure_name", "theory", "heavy_atom_count", "rmsd_angstrom"])
    expected = {(s, t) for s in refs for t in METHODS if t != REFERENCE}
    require(len(rmsd) == len(expected) and {(r["structure_name"], r["theory"]) for r in rmsd} == expected,
            "RMSD table has missing, duplicate, or unknown structure/method pairs")
    for row in rmsd:
        atoms = refs[row["structure_name"]]
        require(int(row["heavy_atom_count"]) == sum(a[0] != "H" for a in atoms), f"Wrong heavy-atom count: {row}")
        require(number(row["rmsd_angstrom"]) >= 0, f"Negative RMSD: {row}")

    bonds = read_table(root / "results/geometry/ts_bond_distance.csv",
                       ["structure_name", "atom_i_1based", "atom_j_1based", "atom_i_symbol",
                        "atom_j_symbol", "distance_angstrom", "theory"])
    seen = set()
    bond_sets = {t: set() for t in METHODS}
    for row in bonds:
        s, theory = row["structure_name"], row["theory"]
        require(s in refs, f"Unknown bond structure: {s}")
        key = (s, theory, int(row["atom_i_1based"]), int(row["atom_j_1based"]))
        require(key not in seen, f"Duplicate bond: {row}")
        seen.add(key)
        require(theory in METHODS, f"Unknown bond theory: {theory}")
        check_bond(row, geometries[METHODS[theory]][s])
        bond_sets[theory].add((s, key[2], key[3]))
    check_bond_coverage(bond_sets, set(refs))
    check_links(root)
    print(f"PASS: {len(cases)} cases, {len(xyz)} XYZ files, {sum(map(len, xyz.values()))} frames, "
          f"{len(folders)} method sets, {len(rmsd)} RMSD rows, {len(bonds)} verified bond rows; local links valid.")
    print("NOTE: ja1g_major is a rotational TS without an assigned forming/breaking bond; see results/geometry/README.md.")


if __name__ == "__main__":
    try:
        validate()
    except (ValueError, OSError, KeyError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        sys.exit(1)
