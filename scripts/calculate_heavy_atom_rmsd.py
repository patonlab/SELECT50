#!/usr/bin/env python3
"""Calculate heavy-atom RMSDs for the deposited SELECT50 geometries.

Expected layout: optimized-geometries/<method>/<structure_id>/<structure_id>.xyz
Reference: wB97X-D3(BJ)-def2-TZVP. Atoms are paired by XYZ row order, hydrogen
is excluded, and all remaining atoms enter an unweighted Kabsch fit without
outlier rejection. Requires NumPy; does not call PyMOL.

Run from any working directory:
    python3 /path/to/SELECT50/scripts/calculate_heavy_atom_rmsd.py --strict

The default output is results/geometry/heavy_atom_rmsd_recalculated.csv, with
source paths and per-comparison status. It leaves heavy_atom_rmsd.csv untouched.
Use --output to choose a different CSV location.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class XYZGeometry:
    symbols: Tuple[str, ...]
    coordinates: np.ndarray


OUTPUT_FIELDS = (
    "structure_name",
    "theory",
    "reference_theory",
    "heavy_atom_count",
    "rmsd_angstrom",
    "theory_xyz",
    "reference_xyz",
    "status",
    "error",
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFERENCE = "wB97X-D3(BJ)-def2-TZVP"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate heavy-atom Kabsch-aligned RMSDs for theory geometries "
            "relative to a reference theory."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--root",
        type=Path,
        default=REPOSITORY_ROOT / "optimized-geometries",
        help="Directory containing the deposited method folders.",
    )
    parser.add_argument(
        "--reference",
        default=DEFAULT_REFERENCE,
        help="Reference theory directory name.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPOSITORY_ROOT / "results/geometry/heavy_atom_rmsd_recalculated.csv",
        help="Output CSV.",
    )
    parser.add_argument(
        "--theories",
        nargs="*",
        default=None,
        help=(
            "Optional theory directories to compare. If omitted, all "
            "directories under --root except the reference are used."
        ),
    )
    parser.add_argument(
        "--include-reference",
        action="store_true",
        help="Also write the reference itself with RMSD = 0.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Return nonzero exit status if any comparison fails.",
    )

    return parser.parse_args()


def normalize_element(symbol: str) -> str:
    """Normalize an XYZ element symbol for comparisons."""
    symbol = symbol.strip()

    # Handle labels like C1 only if encountered, while preserving normal symbols.
    match = re.match(r"^([A-Za-z]{1,2})", symbol)
    if not match:
        return symbol

    element = match.group(1)
    return element[0].upper() + element[1:].lower()


def is_hydrogen(symbol: str) -> bool:
    """Return True for hydrogen (including common H-like labels)."""
    return normalize_element(symbol) == "H"


def read_xyz(path: Path) -> XYZGeometry:
    """Read a standard XYZ/extXYZ geometry."""
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    if not lines:
        raise ValueError("empty XYZ file")

    try:
        natoms = int(lines[0].strip().split()[0])
    except (ValueError, IndexError) as exc:
        raise ValueError("first line does not contain a valid atom count") from exc

    if natoms <= 0:
        raise ValueError(f"invalid atom count: {natoms}")

    if len(lines) < natoms + 2:
        raise ValueError(
            f"XYZ header reports {natoms} atoms, but file is too short"
        )

    if any(line.strip() for line in lines[natoms + 2:]):
        raise ValueError("expected a single XYZ geometry; extra frames or trailing data found")

    symbols: List[str] = []
    coords: List[List[float]] = []

    for atom_index, line in enumerate(lines[2 : 2 + natoms], start=1):
        fields = line.split()

        if len(fields) < 4:
            raise ValueError(
                f"atom line {atom_index} has fewer than four fields"
            )

        symbol = normalize_element(fields[0])

        try:
            xyz = [float(fields[1]), float(fields[2]), float(fields[3])]
        except ValueError as exc:
            raise ValueError(
                f"invalid coordinates on atom line {atom_index}"
            ) from exc

        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(
                f"non-finite coordinates on atom line {atom_index}"
            )

        symbols.append(symbol)
        coords.append(xyz)

    return XYZGeometry(
        symbols=tuple(symbols),
        coordinates=np.asarray(coords, dtype=float),
    )


def validate_atom_mapping(
    reference: XYZGeometry,
    mobile: XYZGeometry,
) -> None:
    """
    Require identical atom count and atom ordering.

    This is intentionally strict: an RMSD is only meaningful here if atom i in
    one geometry corresponds to atom i in the other geometry.
    """
    if len(reference.symbols) != len(mobile.symbols):
        raise ValueError(
            f"atom-count mismatch: reference={len(reference.symbols)}, "
            f"theory={len(mobile.symbols)}"
        )

    mismatches = [
        (i + 1, ref, mob)
        for i, (ref, mob) in enumerate(
            zip(reference.symbols, mobile.symbols)
        )
        if ref.casefold() != mob.casefold()
    ]

    if mismatches:
        preview = ", ".join(
            f"{i}:{ref}/{mob}"
            for i, ref, mob in mismatches[:10]
        )

        suffix = (
            ""
            if len(mismatches) <= 10
            else f" ... ({len(mismatches)} mismatches total)"
        )

        raise ValueError(
            "atom-order/element mismatch at "
            + preview
            + suffix
        )


def heavy_atom_coordinates(
    reference: XYZGeometry,
    mobile: XYZGeometry,
) -> Tuple[np.ndarray, np.ndarray, int]:
    validate_atom_mapping(reference, mobile)

    heavy_mask = np.asarray(
        [not is_hydrogen(symbol) for symbol in reference.symbols],
        dtype=bool,
    )

    n_heavy = int(np.sum(heavy_mask))

    if n_heavy == 0:
        raise ValueError("structure contains no heavy atoms")

    return (
        reference.coordinates[heavy_mask].copy(),
        mobile.coordinates[heavy_mask].copy(),
        n_heavy,
    )


def kabsch_rmsd(
    reference_coords: np.ndarray,
    mobile_coords: np.ndarray,
) -> float:
    """
    Return RMSD after optimal rigid-body alignment with the Kabsch algorithm.

    Both arrays must have shape (N, 3), with corresponding atoms in the same
    order.
    """
    if reference_coords.shape != mobile_coords.shape:
        raise ValueError(
            f"coordinate-shape mismatch: {reference_coords.shape} vs "
            f"{mobile_coords.shape}"
        )

    if (
        reference_coords.ndim != 2
        or reference_coords.shape[1] != 3
    ):
        raise ValueError(
            f"coordinates must have shape (N, 3), got {reference_coords.shape}"
        )

    if len(reference_coords) == 0:
        raise ValueError("no coordinates supplied")

    if not np.all(np.isfinite(reference_coords)) or not np.all(np.isfinite(mobile_coords)):
        raise ValueError("coordinates must be finite")

    # Translate both structures to their heavy-atom centroids.
    ref_centered = (
        reference_coords
        - reference_coords.mean(axis=0)
    )
    mob_centered = (
        mobile_coords
        - mobile_coords.mean(axis=0)
    )

    # Find the optimal rotation that maps mobile onto reference.
    covariance = mob_centered.T @ ref_centered
    U, _singular_values, Vt = np.linalg.svd(covariance)

    rotation = U @ Vt

    # Prevent an improper rotation/reflection.
    if np.linalg.det(rotation) < 0:
        U[:, -1] *= -1.0
        rotation = U @ Vt

    aligned_mobile = mob_centered @ rotation

    difference = aligned_mobile - ref_centered

    rmsd = np.sqrt(
        np.sum(difference * difference)
        / len(reference_coords)
    )

    return float(rmsd)


def choose_xyz(structure_dir: Path, structure_name: str) -> Optional[Path]:
    """Use the exact current filename; never substitute an unrelated geometry."""
    path = structure_dir / f"{structure_name}.xyz"
    return path if path.is_file() else None


def find_structure_directories(
    theory_root: Path,
) -> dict[str, Path]:
    """
    Return direct child structure directories keyed case-insensitively.

    The benchmark directory layout used here places one structure directory
    directly below each theory directory.
    """
    result: dict[str, Path] = {}

    if not theory_root.is_dir():
        return result

    for path in theory_root.iterdir():
        if not path.is_dir():
            continue

        name_key = path.name.casefold()

        result[name_key] = path

    return result


def write_row(
    writer: csv.DictWriter,
    structure_name: str,
    theory: str,
    reference: str,
    heavy_atom_count: str = "",
    rmsd: str = "",
    theory_xyz: str = "",
    reference_xyz: str = "",
    status: str = "",
    error: str = "",
) -> None:
    writer.writerow(
        {
            "structure_name": structure_name,
            "theory": theory,
            "reference_theory": reference,
            "heavy_atom_count": heavy_atom_count,
            "rmsd_angstrom": rmsd,
            "theory_xyz": theory_xyz,
            "reference_xyz": reference_xyz,
            "status": status,
            "error": error,
        }
    )


def main() -> int:
    args = parse_args()

    root = args.root.expanduser().resolve()
    output = args.output.expanduser().resolve()

    if not root.is_dir():
        print(
            f"ERROR: root directory does not exist: {root}",
            file=sys.stderr,
        )
        return 2

    reference_root = root / args.reference

    if not reference_root.is_dir():
        print(
            f"ERROR: reference directory does not exist: "
            f"{reference_root}",
            file=sys.stderr,
        )
        return 2

    if args.theories:
        theories = list(args.theories)
    else:
        theories = sorted(
            [
                path.name
                for path in root.iterdir()
                if path.is_dir()
                and path.name.casefold()
                != args.reference.casefold()
            ],
            key=str.casefold,
        )

    reference_dirs = find_structure_directories(reference_root)

    if not reference_dirs:
        print(
            f"ERROR: no structure directories found below "
            f"{reference_root}",
            file=sys.stderr,
        )
        return 2

    theory_dirs = {
        theory: find_structure_directories(root / theory)
        for theory in theories
    }

    rows_written = 0
    success_count = 0
    failure_count = 0

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output.open(
        "w",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_FIELDS,
        )
        writer.writeheader()

        for ref_key, ref_structure_dir in sorted(
            reference_dirs.items(),
            key=lambda item: item[1].name.casefold(),
        ):
            structure_name = ref_structure_dir.name

            reference_xyz = choose_xyz(
                ref_structure_dir,
                structure_name,
            )

            if reference_xyz is None:
                message = "no usable reference XYZ found"

                for theory in theories:
                    write_row(
                        writer,
                        structure_name,
                        theory,
                        args.reference,
                        status="error",
                        error=message,
                    )
                    rows_written += 1
                    failure_count += 1

                continue

            try:
                reference_geometry = read_xyz(reference_xyz)
            except (OSError, ValueError) as exc:
                message = f"reference XYZ error: {exc}"

                for theory in theories:
                    write_row(
                        writer,
                        structure_name,
                        theory,
                        args.reference,
                        reference_xyz=str(reference_xyz),
                        status="error",
                        error=message,
                    )
                    rows_written += 1
                    failure_count += 1

                continue

            if args.include_reference:
                heavy_count = sum(
                    not is_hydrogen(symbol)
                    for symbol in reference_geometry.symbols
                )

                write_row(
                    writer,
                    structure_name,
                    args.reference,
                    args.reference,
                    heavy_atom_count=str(heavy_count),
                    rmsd="0.000000",
                    theory_xyz=str(reference_xyz),
                    reference_xyz=str(reference_xyz),
                    status="ok",
                )
                rows_written += 1
                success_count += 1

            for theory in theories:
                theory_root = root / theory

                if not theory_root.is_dir():
                    write_row(
                        writer,
                        structure_name,
                        theory,
                        args.reference,
                        reference_xyz=str(reference_xyz),
                        status="error",
                        error=f"theory directory not found: {theory_root}",
                    )
                    rows_written += 1
                    failure_count += 1
                    continue

                structure_dir = theory_dirs[theory].get(ref_key)

                if structure_dir is None:
                    write_row(
                        writer,
                        structure_name,
                        theory,
                        args.reference,
                        reference_xyz=str(reference_xyz),
                        status="error",
                        error="matching structure directory not found",
                    )
                    rows_written += 1
                    failure_count += 1
                    continue

                theory_xyz = choose_xyz(
                    structure_dir,
                    structure_name,
                )

                if theory_xyz is None:
                    write_row(
                        writer,
                        structure_name,
                        theory,
                        args.reference,
                        reference_xyz=str(reference_xyz),
                        status="error",
                        error="no usable theory XYZ found",
                    )
                    rows_written += 1
                    failure_count += 1
                    continue

                try:
                    theory_geometry = read_xyz(theory_xyz)

                    (
                        reference_heavy,
                        theory_heavy,
                        heavy_count,
                    ) = heavy_atom_coordinates(
                        reference_geometry,
                        theory_geometry,
                    )

                    rmsd = kabsch_rmsd(
                        reference_heavy,
                        theory_heavy,
                    )

                    write_row(
                        writer,
                        structure_name,
                        theory,
                        args.reference,
                        heavy_atom_count=str(heavy_count),
                        rmsd=f"{rmsd:.6f}",
                        theory_xyz=str(theory_xyz),
                        reference_xyz=str(reference_xyz),
                        status="ok",
                    )

                    success_count += 1

                except (OSError, ValueError, np.linalg.LinAlgError) as exc:
                    write_row(
                        writer,
                        structure_name,
                        theory,
                        args.reference,
                        theory_xyz=str(theory_xyz),
                        reference_xyz=str(reference_xyz),
                        status="error",
                        error=str(exc),
                    )
                    failure_count += 1

                rows_written += 1

    print(f"Root:                    {root}")
    print(f"Reference theory:        {args.reference}")
    print(f"Reference structures:    {len(reference_dirs)}")
    print(f"Theories compared:       {len(theories)}")
    print(f"Successful RMSDs:        {success_count}")
    print(f"Failed comparisons:      {failure_count}")
    print(f"Rows written:            {rows_written}")
    print(f"Output CSV:              {output}")

    if theories:
        print("\nTheories:")
        for theory in theories:
            print(f"  {theory}")

    if args.strict and failure_count:
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
