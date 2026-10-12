# Palladium-cyclization landscape — Figure 7

This study maps to manuscript Figure 7 and SI Table S21/Figure S14 in the separately supplied working drafts. The literature source is printed as main reference 82 and SI reference 80, DOI `10.1021/acscatal.3c00732`. This guide describes deposited files; it does not establish missing calculation provenance or resolve manuscript discrepancies.

## Files and units

| Material | Location and contents |
| --- | --- |
| UMA stationary structures | [structures/uma/](structures/uma/): 21 individual XYZs and [all_21_frames.xyz](structures/uma/all_21_frames.xyz) |
| Electronic/Gibbs energies | [gibbs_21_frames.csv](gibbs_21_frames.csv): 21 species, absolute E/G and G−E in hartree; recorded temperature 298.15 K and pressure 1 atm; imaginary frequencies in cm⁻¹ |
| Spectator-adjusted profile | [profile_species.csv](profile_species.csv): absolute `E`, `G`, `E_total`, `G_total` in hartree; `G_rel_kcal` in kcal/mol relative to the separated reference row `3a`; `spectators_added` specifies coefficients |
| IRC paths | [scan_paths/irc/](scan_paths/irc/): TS1, TS1′, TS2, TS3, TS3′, TS4, TS4′ coordinate/energy records |
| Relaxed dihedral scan | [scan_dihedral.csv](scan_paths/scan_dihedral/scan_dihedral.csv): 73 points, degrees, absolute E in hartree and relative E in kcal/mol; frame paths and discontinuity flags |
| NEB/MEP paths | [scan_paths/neb/](scan_paths/neb/): TS1, TS4, TS4′ bands; absolute E in hartree, relative E in kcal/mol, path distances in Å where recorded |
| Published DFT structures | [structures/published_dft/](structures/published_dft/): 22 XYZs, [manifest](structures/published_dft/manifest.csv), and [extraction notes](structures/published_dft/README.txt) |
| Geometry comparisons | [table_s20_equivalent_pd.csv](table_s20_equivalent_pd.csv), [pd_distances.csv](pd_distances.csv), [pd_angles.csv](pd_angles.csv): RMSD/distances in Å, angles in degrees, one-based atom indices |

XYZ coordinates use Å. Helper/spectator rows in the profile have blank adjusted/relative energies where a complete profile state is not defined. The energy CSV with explicit hartree headers rounds some absolute values more than `profile_species.csv`; retain the distinction when recreating a profile.

## Figure and label mapping

- Figure 7A uses spectator-adjusted relative Gibbs energies. `3a` corresponds to SI `Sub`; `4a` to `Prod`; `4a_prime` is the second product. Folder/CSV `Int1`, `int2`–`int5`, and `ts1`–`ts4` correspond to the formatted SI labels; `_prime` represents a prime mark.
- Figure 7B/S14 use IRC paths. The requested unprimed TS1–TS4 and TS4′ energy tables contain 37, 16, 11, 51, and 51 points, respectively. File presence alone does not verify relaxed endpoint connectivity.
- Figure 7C uses the 73-point UMA dihedral scan and the deposited DFT/UMA overlay table. Raw scan angles and plotting angles are distinct columns. The source/digitization procedure for the published DFT scan and treatment of discontinuities need documentation.
- Figure 7D uses product-forming NEB paths. Deposited bands contain 18 frames for TS1 and 22 each for TS4/TS4′. The SI's general statement of 16 intermediate images must be reconciled per path, distinguishing endpoints and separately sampled MEP points.
- Figure 7E compares four TS overlays. The table reproduces the quoted all-atom RMSDs: TS1 1.031, TS3 0.446, TS4 0.353, TS4′ 0.513 Å. Heavy-atom RMSDs are separate columns.

## Outstanding provenance and document corrections

The working SI describes ORCA 6.1.1, UMA-S-1.1/OMol through an external-method interface, RTX 5070, numerical frequencies, TIGHTOPT intermediate optimization, and bidirectional IRCs. These are reported settings; original per-job inputs/model records remain needed to verify each deposited calculation and thermochemistry/standard-state treatment.

S21 currently has 20 species and copies product E/G into Int4. The deposited Int4 values are E = −3522.466642034 and G = −3521.58878451 hartree, whereas the SI lists −1015.926108778 and −1015.63207358. `4a_prime` is also absent from S21. The documents have not been corrected during the guide/metadata work.

The Pd contact comparison tables select contacts within 2.8 Å. They are not the exact set of distances/angles displayed in Figure 7E, which also includes forming bonds. Deposit those specific atom selections and fit/atom-correspondence conventions before treating the broad table as proof of every selected-metric claim.

DFT coordinates are reported as extracted from Section 26 of `cs3c00732_si_001.pdf`; the manifest records page starts and Eopt in hartree. The PDF is not deposited. Original structure identification, extraction/atom mapping, and source notices need a complete provenance account.

The main draft compares the calculation with related experimental substrates at 60 °C, explicitly noting that the exact modeled system was not experimentally characterized. Per-entry ee values, source pages, and the conversion account are still missing. The calculated 298.15 K Gibbs energies and experimental selectivity at 333.15 K are different conditions.
