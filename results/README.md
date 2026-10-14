# SELECT50 benchmark data

The core dataset contains 50 transition-state pairs, with reference structures, method-specific geometries, relative energies and structural comparisons. Reaction landscapes, large systems and experimental calibration are separate datasets.

## Files

| Location | Contents |
| --- | --- |
| [reference/starting_structures/](reference/starting_structures/) | 100 reference XYZ files, identical to the corresponding optimized DFT geometries |
| [reference/](reference/README.md) | Case metadata, energy orientations and calculation exclusions |
| [optimized-geometries/](../optimized-geometries/) | Nine method sets: seven with 100 structures each, and g-xTB/r²SCAN-3c with 98 each |
| [energetics/](energetics/) | Two 50-case master tables, raw electronic/Gibbs energies and energy summaries |
| [geometry/](geometry/README.md) | 796 heavy-atom RMSDs, 1,488 bond distances, 1,322 bond errors and summary tables |
| [perturbed_guesses/ts_recovery_full_data.csv](perturbed_guesses/ts_recovery_full_data.csv) | 300 recovery runs; calculation records are in [ts_recovery/](../ts_recovery/) |

## Energies

[delta_delta_e.csv](energetics/delta_delta_e.csv) and [delta_delta_g.csv](energetics/delta_delta_g.csv) contain signed minor-minus-major gaps in kcal/mol. Major is the lower reference electronic-energy structure; the same role is used across methods and Gibbs energies. Five reference Gibbs gaps are negative. g-xTB evaluates 49 pairs, with an82 blank in both master tables.

Six `*_raw_energies.csv` files contain 100 electronic energies each for DFT, eSEN, g-xTB, MACE-OMOL, UMA-M and UMA-S. Five `*_raw_Gibbs_energies.csv` files contain Gibbs energies for DFT, eSEN, g-xTB, UMA-M and UMA-S. MACE electronic energies, thermal corrections and Gibbs energies are recorded together in [mace_omol_raw_thermochemistry.csv](energetics/mace_omol_raw_thermochemistry.csv). Absolute energies are in hartree and structure IDs use `<case_id>_major` or `<case_id>_minor`.

The energy summary CSVs report overall MAE/RMSE, category MAEs and ordering accuracy. The [validator](../scripts/validate_data.py) checks these against the master tables and checks supported raw gaps within 0.001 kcal/mol, accommodating the precision of the absolute-energy records.

## Geometry

Heavy-atom RMSDs use fixed XYZ atom order and an equally weighted proper Kabsch rotation. Bond indices are one-based and distances use ångströms. Both ja1g structures are excluded from bond analysis, giving 166 reference bonds across 98 structures. Both an82 structures are excluded from g-xTB and r²SCAN-3c geometry comparisons. See the [geometry guide](geometry/README.md).

## Related datasets and figure inputs

| Dataset | Main figure | SI tables | SI figures |
| --- | --- | --- | --- |
| [Cycloaddition landscape](../reaction_landscapes/cycloaddition/) | 5 | S19–S20 | S13 |
| [Palladium-cyclization landscape](../reaction_landscapes/palladium_cyclization/) | 6 | S21 | S14 |
| [Large systems](../large_systems/) | 7 | S22 | — |
| [Ni calibration](../calibration/nickel_epoxide_arylation/) | 8A | S23–S24 | S15 |
| [Pd calibration](../calibration/palladium_enone_addition/) | 8B | S23–S24 | S16 |

The recovery calculations use UMA-S-1.2 with the OMol task. The two landscape datasets use UMA-S-1.1 with the OMol task. Large-system and calibration calculations use UMA-S-1.2. Ni/Pd experimental gaps use 307.15 K and 333 K, respectively.
