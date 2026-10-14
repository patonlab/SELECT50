# Geometry comparisons

## Files and units

| File | Contents |
| --- | --- |
| [heavy_atom_rmsd.csv](heavy_atom_rmsd.csv) | 796 heavy-atom RMSDs relative to the DFT geometry reference |
| [ts_bond_distances.csv](ts_bond_distances.csv) | 1,488 distances with one-based atom indices |
| [ts_bond_errors.csv](ts_bond_errors.csv) | 1,322 signed and absolute non-reference bond errors |
| [heavy_atom_rmsd_method_summary.csv](heavy_atom_rmsd_method_summary.csv) | Mean and median RMSD by method |
| [ts_bond_method_summary.csv](ts_bond_method_summary.csv) | Bond RMSE and MAE by method |
| [ts_bond_category_summary.csv](ts_bond_category_summary.csv) | Bond errors grouped by chemical category |

Distances and geometric errors use ångströms. Both ja1g structures are excluded from bond analysis, giving 166 reference bonds across 98 structures. Each of the six complete comparison methods has 166 bonds; g-xTB and r²SCAN-3c each have 163 bonds across 96 structures because both an82 geometries are excluded. Both ja1g structures remain in the RMSD dataset.

## RMSD definition

The [RMSD script](../../scripts/calculate_heavy_atom_rmsd.py) pairs computed and reference atoms by XYZ row order, verifies atom counts and element sequences, and excludes hydrogen. It aligns all remaining atoms with equal weight using the Kabsch algorithm and a proper rotation:

$$
\mathrm{RMSD}=\min_{\mathbf{R},\mathbf{t}}
\sqrt{\frac{1}{N}\sum_{i=1}^{N}
\left\|\mathbf{x}_i\mathbf{R}+\mathbf{t}-\mathbf{y}_i\right\|^2},
\qquad \mathbf{R}^{\mathsf T}\mathbf{R}=\mathbf{I},\quad \det(\mathbf{R})=+1.
$$

Here, `N` is the number of heavy atoms. The fit uses centroid subtraction and singular value decomposition, with no mass weighting, outlier rejection or symmetry-based atom remapping. The reference geometries are ωB97X-D3(BJ)/def2-TZVP. Separate landscape and large-system tables report all-atom RMSDs where indicated.

## Recalculate

From the repository root:

```sh
python3 -m pip install numpy
python3 scripts/calculate_heavy_atom_rmsd.py --output rmsd_recalculated.csv
```

The script compares each structure across the eight non-reference methods. Its output includes source paths, heavy-atom counts, RMSDs and calculation status. `--output` selects the destination; the default is `results/geometry/heavy_atom_rmsd_recalculated.csv`. Existing output at the selected path is overwritten. Use `-h` for all options; `--strict` returns a nonzero exit code if any comparison fails.
