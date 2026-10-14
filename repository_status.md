# Repository status

Updated 2026-10-14 to describe the current deposited files and analysis conventions.

## Core inventory

- 50 cases including an82; 100 reference copies, each byte-identical to the optimized DFT reference.
- 896 optimized XYZs: seven methods ×100 and g-xTB/r²SCAN-3c ×98.
- Two 50-row pair-energy tables; g-xTB evaluates 49 pairs in each.
- Raw core energies in results/energetics/: six 100-row electronic-energy tables, five Gibbs-energy tables and the MACE-OMOL thermochemistry table.
- 796 heavy-atom RMSDs, 1,488 bond distances and 1,322 non-reference bond errors.
- 166 reference bonds over 98 structures. Both ja1g structures are excluded from bond analysis. g-xTB/r²SCAN-3c each cover 163 bonds over 96 structures.
- 600 source-orientation records, case metadata and explicit geometry-failure ledgers.

Major/minor identities are assigned using the reference DFT electronic energies and retained for Gibbs comparisons. Five reference Gibbs gaps are negative. Core model labels are UMA-M-1.1 and UMA-S-1.2. MACE electronic ordering is 46/50 after the ac27, d17c and ja03 sign corrections.

The ac16 and ac30 DOI fields are synchronized across case metadata and both pair-energy tables: `10.1021/acscatal.1c00616` and `10.1021/acsomega.9b00030`, respectively.

## Validation

Run `python scripts/validate_data.py`. It checks XYZ parsing, method/case coverage, reference-copy equality, element ordering, metadata atom counts, CSV schemas, declared energy exclusions, E/G orientation records, energy method summaries, bond distances/errors and local Markdown links.

Tests are available with `python -m unittest discover -s tests`; the RMSD tests require NumPy. The [RMSD calculator](scripts/calculate_heavy_atom_rmsd.py) uses fixed XYZ atom order, heavy atoms and a proper Kabsch rotation. Validation of deposited tables is distinct from authentication of original calculation outputs.

## Related datasets and figure map

| Dataset | Main figure | SI tables/figures | Deposited records |
| --- | --- | --- | --- |
| Cycloaddition landscape | 5 | S19/S20; S13 | Stationary structures, E/G, 2D scan, IRC and activation-strain tables |
| Palladium-cyclization landscape | 6 | S21; S14 | 20 UMA stationary structures, 21 published DFT structures, profiles, IRC, dihedral scan and NEB/MEP paths |
| Large systems | 7 | S22 | 12 optimized/start pairs and an E/G/frequency/timing table |
| Ni/Pd calibration | 8A/B | S23/S24; S15/S16 | 13/15 entries, 100 partitions each, memberships, coefficients and held-out predictions |

The recovery dataset contains 300 logs and 300 XYZ exports. Its [300-row summary](results/perturbed_guesses/ts_recovery_full_data.csv) records recovery metrics and classifications. The recovery XYZ count includes the starting-geometry export for jo8v_min_rmsd0p14_r05.

Calibration experimental gaps use the author-confirmed temperatures: Ni 307.15 K and Pd 333 K. The refit retains the existing 100 partitions, seed record and raw UMA predictions. Partition plotting files provide raw/calibrated MAEs and separately labeled mean-prediction null errors; per-entry plotting files provide Figure 8 coordinates and sample SDs. The shared [refresh driver](calibration/refresh_saved_calibration.py) and family-level `calibration_refit.json` files record this conversion and refit.
