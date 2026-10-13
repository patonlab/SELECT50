# Core SELECT50 benchmark data

The core contains [100 reference XYZ copies](reference/starting_structures/), two [50-case energy tables](energetics/), and [geometry comparisons](geometry/). The nine [optimized geometry sets](../optimized-geometries/) contain 896 structures: seven sets of 100 and g-xTB/r²SCAN-3c sets of 98. The latter exclude both failed an82 structures. The obsolete ja30 copies are retained only in the [historical archive](../archive/replaced_ja30/README.md).

Major is defined by lower reference electronic energy and retained for both E and G. Energy differences are signed minor minus major in kcal/mol. Five reference G differences are negative. g-xTB has 49 evaluated energy pairs and an explicitly blank failed an82 value in each table. [Metadata and orientations](reference/README.md) preserve the source-to-current mapping.

There are 796 heavy-atom RMSDs, 1,488 bond distances and 1,322 non-reference bond errors. Both ja1g structures are excluded from bond analysis by the author's 2026-10-12 confirmation: 166 reference bonds cover 98 structures; g-xTB/r²SCAN-3c each have 163 bonds over 96 structures. Atom indices are one-based and distances use ångströms.

The [RMSD script](../scripts/calculate_heavy_atom_rmsd.py) uses fixed XYZ atom order, excludes hydrogen, and applies an equally weighted proper Kabsch rotation without outlier rejection. All 796 available results reproduce at six-decimal output precision; four absent an82 comparisons are reported explicitly. See [geometry methods](geometry/README.md).

Recovery logs and case-study/calibration data are deposited elsewhere in this repository. Original reference frequency/single-point outputs and complete per-job checkpoint/environment provenance remain incomplete. [Repository status](../repository_status.md) separates current validation from those outstanding requirements.
