# Geometry comparisons

[heavy_atom_rmsd.csv](heavy_atom_rmsd.csv) contains 796 heavy-atom RMSDs in ångströms. [ts_bond_distances.csv](ts_bond_distances.csv) contains 1,488 distances: 166 assigned reference bonds across 98 structures, with both ja1g structures excluded. The six complete comparison methods have 166 bonds each; g-xTB and r²SCAN-3c have 163 each over 96 structures because both an82 calculations failed. [ts_bond_errors.csv](ts_bond_errors.csv) contains the 1,322 non-reference comparisons. Atom indices are one-based.

## RMSD methodology

The [analysis script](../../scripts/calculate_heavy_atom_rmsd.py) calculates the core heavy-atom RMSDs directly with NumPy and the Kabsch algorithm.

### Atom correspondence and selection

Each computed structure is compared with its identically named reference under `optimized-geometries/wB97X-D3(BJ)-def2-TZVP/`. Atom correspondence is fixed by XYZ row order. The script rejects unequal atom counts or differing element sequences, then excludes all hydrogen atoms from both geometries. Every remaining atom, including metals, enters the fit with equal weight. No atoms are removed as outliers and no symmetry-equivalent atom permutations are searched. Matching element sequences verifies consistency but does not by itself prove chemical correspondence if atoms of the same element were reordered.

### Mathematical definition and Kabsch fit

For `N` heavy atoms, let $\mathbf{x}_i$ and $\mathbf{y}_i$ be row vectors containing the computed and reference coordinates. The fitted RMSD is

$$
\mathrm{RMSD}=\min_{\mathbf{R},\mathbf{t}}
\sqrt{\frac{1}{N}\sum_{i=1}^{N}
\left\|\mathbf{x}_i\mathbf{R}+\mathbf{t}-\mathbf{y}_i\right\|^2},
\qquad \mathbf{R}^{\mathsf T}\mathbf{R}=\mathbf{I},\quad \det(\mathbf{R})=+1.
$$

The arithmetic centroids are $\bar{\mathbf{x}}=N^{-1}\sum_i\mathbf{x}_i$ and $\bar{\mathbf{y}}=N^{-1}\sum_i\mathbf{y}_i$. The script subtracts these centroids to form the centered coordinate matrices $\mathbf{X}_c$ and $\mathbf{Y}_c$. It computes the covariance $\mathbf{H}=\mathbf{X}_c^{\mathsf T}\mathbf{Y}_c$ and its singular value decomposition $\mathbf{H}=\mathbf{U}\mathbf{\Sigma}\mathbf{V}^{\mathsf T}$. The optimal proper rotation is

$$
\mathbf{R}=\mathbf{U}\mathbf{D}\mathbf{V}^{\mathsf T},\qquad
\mathbf{D}=\operatorname{diag}\!\left(1,1,\operatorname{sign}\det(\mathbf{U}\mathbf{V}^{\mathsf T})\right).
$$

The determinant correction prevents reflection. The script implements it by reversing the final column of $\mathbf{U}$ when necessary. After rotation, it calculates the RMSD from $\mathbf{X}_c\mathbf{R}-\mathbf{Y}_c$; the corresponding translation for the original coordinates is $\mathbf{t}=\bar{\mathbf{y}}-\bar{\mathbf{x}}\mathbf{R}$.

This is a rigid fit: internal bond lengths are unchanged and no geometry optimization occurs. Each heavy atom contributes equally, so the RMSD is not mass weighted. The denominator is the atom count `N`, not `3N`. Results are expressed in ångströms.

### Recalculate the benchmark

Python and NumPy are required. From the repository root:

```sh
python3 -m pip install numpy
python3 scripts/calculate_heavy_atom_rmsd.py --output /tmp/select50_rmsd.csv
```

The script resolves its default paths relative to its own location, so it also works when launched from another directory. It uses the current `<method>/<structure_id>/<structure_id>.xyz` layout, including `MACE-OMOL`, and compares all 100 structures across the eight non-reference methods.

By default it writes `results/geometry/heavy_atom_rmsd_recalculated.csv`, leaving the deposited `heavy_atom_rmsd.csv` unchanged. The recalculated CSV includes method/reference folder names, heavy-atom counts, RMSDs, source paths, and success/error status. Rerunning overwrites that recalculated file. `--output` selects another destination; `--strict` returns a nonzero exit status if any comparison fails. Missing or invalid geometries are reported rather than replaced by another XYZ file.

### Relationship to PyMOL `align` and `cycles=0`

In PyMOL, `align` with `cycles=0` disables iterative outlier rejection while preserving superposition; a heavy-atom selection must separately exclude hydrogen. PyMOL performs its own atom-matching procedure, so agreement with a fixed-index Kabsch fit requires the same atom pairs and selections. See the [PyMOL command reference](https://pymol.org/pymol-command-ref.html#align).

`cycles=0` is not a parameter of the deposited NumPy script. Its equivalent methodological property is that every selected heavy atom enters a single Kabsch fit without outlier rejection. Supporting tables that explicitly report all-atom RMSDs include hydrogen and must retain that distinction.

### Methods text for the core heavy-atom comparisons

> Heavy-atom root-mean-square deviations (RMSDs) were calculated relative to the corresponding ωB97X-D3(BJ)/def2-TZVP reference geometries using the deposited NumPy analysis script. Atoms were paired by their order in the XYZ files after verifying identical atom counts and element sequences. Hydrogen atoms were excluded, and all remaining atoms were retained with equal weight. Both coordinate sets were centered on their heavy-atom centroids and subjected to Kabsch least-squares superposition. A determinant correction restricted the transformation to a proper rotation. No outlier rejection or symmetry-based atom remapping was applied. RMSDs were calculated as the square root of the mean squared Cartesian separation after superposition and are reported in ångströms.

## Bond-analysis definitions

### ac32_minor

The r2SCAN-3c C17–H25 distance is **2.354753 Å** in [the deposited geometry](../../optimized-geometries/r2SCAN-3c/ac32_minor/ac32_minor.xyz).

### ja1g

Both `ja1g_major` and `ja1g_minor` are excluded from bond analysis. The bond dataset contains 166 reference bonds across 98 structures. Both ja1g structures are included in the coordinate and RMSD datasets. See [gallery sheet 2](../../assets/figure_02.png) and [roles and exclusions](../reference/README.md#roles-and-exclusions).
