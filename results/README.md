# Core SELECT50 benchmark data

The core contains [100 reference XYZ copies](reference/starting_structures/), two [50-case energy tables](energetics/), and [geometry comparisons](geometry/). The nine [optimized geometry sets](../optimized-geometries/) contain 896 structures: seven sets of 100 and g-xTB/r²SCAN-3c sets of 98. The latter exclude both failed an82 structures.

Major is defined by lower reference electronic energy and retained for both E and G. Energy differences are signed minor minus major in kcal/mol. Five reference G differences are negative. g-xTB has 49 evaluated energy pairs and an explicitly blank failed an82 value in each table. [Metadata and orientations](reference/README.md) preserve the source-to-current mapping.

There are 796 heavy-atom RMSDs, 1,488 bond distances and 1,322 non-reference bond errors. Both ja1g structures are excluded from bond analysis: 166 reference bonds cover 98 structures; g-xTB/r²SCAN-3c each have 163 bonds over 96 structures. Atom indices are one-based and distances use ångströms.

The [RMSD script](../scripts/calculate_heavy_atom_rmsd.py) uses fixed XYZ atom order, excludes hydrogen, and applies an equally weighted proper Kabsch rotation without outlier rejection. See [geometry methods](geometry/README.md).

The [recovery summary](perturbed_guesses/ts_recovery_full_data.csv) contains 300 runs; the corresponding logs and XYZ exports are in [ts_recovery/](../ts_recovery/). [Reaction landscapes](../reaction_landscapes/), [large-system pairs](../large_systems/large_system_transition_pairs.csv), and [calibration data](../calibration/) are separate datasets.

Six raw-energy tables in [energetics/](energetics/) record 100 electronic energies each for DFT, eSEN, g-xTB, MACE-OMOL, UMA-M and UMA-S. Their `structure_id` fields use the same `<case_id>_major/minor` convention as the core reference files. Raw energies are in hartree; canonical pair differences and summary errors are in kcal/mol. g-xTB's raw an82 entries do not alter its declared benchmark exclusion.
