# Palladium-cyclization landscape — Figure 6

This study maps to manuscript Figure 6 and SI Table S21/Figure S14. The literature source DOI is `10.1021/acscatal.3c00732`. The reaction-landscape records identify UMA-S-1.1 with the OMol task.

## Files and units

| Material | Location and contents |
| --- | --- |
| UMA stationary structures | [structures/uma/](structures/uma/): 20 individual XYZs and [all_20_frames.xyz](structures/uma/all_20_frames.xyz) |
| Electronic/Gibbs energies | [gibbs_20_frames.csv](gibbs_20_frames.csv): 20 species, absolute E/G and G−E in hartree; recorded temperature 298.15 K and pressure 1 atm; imaginary frequencies in cm⁻¹ |
| Spectator-adjusted profile | [profile_species.csv](profile_species.csv): absolute `E`, `G`, `E_total`, `G_total` in hartree; `G_rel_kcal` in kcal/mol relative to the separated reference row `3a`; `spectators_added` specifies coefficients |
| IRC paths | [scan_paths/irc/](scan_paths/irc/): TS1, TS1′, TS2, TS3, TS3′, TS4, TS4′ coordinate/energy records |
| Relaxed dihedral scan | [scan_dihedral.csv](scan_paths/scan_dihedral/scan_dihedral.csv): 73 points, degrees, absolute E in hartree and relative E in kcal/mol; frame paths and discontinuity flags |
| NEB/MEP paths | [scan_paths/neb/](scan_paths/neb/): TS1, TS4, TS4′ bands; absolute E in hartree, relative E in kcal/mol, path distances in Å where recorded |
| Published DFT structures | [structures/published_dft/](structures/published_dft/): 21 XYZs, [manifest](structures/published_dft/manifest.csv), and [extraction notes](structures/published_dft/README.txt) |
| Geometry comparisons | [table_s20_equivalent_pd.csv](table_s20_equivalent_pd.csv), [pd_distances.csv](pd_distances.csv), [pd_angles.csv](pd_angles.csv): RMSD/distances in Å, angles in degrees, one-based atom indices |

XYZ coordinates use Å. Helper/spectator rows in the profile have blank adjusted/relative energies where a complete profile state is not defined. The energy CSV with explicit hartree headers rounds some absolute values more than `profile_species.csv`; retain the distinction when recreating a profile.

## Figure and label mapping

- Figure 6A uses spectator-adjusted relative Gibbs energies. `3a` corresponds to SI `Sub`; `4a` to `Prod`. Folder/CSV `Int1`, `int2`–`int5`, and `ts1`–`ts4` correspond to the formatted SI labels; `_prime` represents a prime mark.
- Figure 6B/S14 use IRC paths. The unprimed TS1–TS4 and TS4′ energy tables contain 37, 16, 11, 51, and 51 points, respectively.
- Figure 6C uses the 73-point UMA dihedral scan and the [DFT/UMA overlay table](scan_paths/scan_dihedral/DFT_UMA_relaxed_scan_overlay.csv). Raw scan angles and plotting angles are distinct columns.
- Figure 6D uses product-forming NEB paths. Deposited bands contain 18 frames for TS1 (16 intermediate images) and 22 each for TS4/TS4′ (20 intermediate images). MEP profiles are separately sampled paths.
- Figure 6E compares four TS overlays. The tabulated all-atom RMSDs round to TS1 1.031, TS3 0.446, TS4 0.353, TS4′ 0.513 Å. Heavy-atom RMSDs are separate columns.

## Structure and contact records

Table S21 lists the 20 retained UMA stationary species. The Pd contact comparison tables select contacts within 2.8 Å and record one-based atom indices for distances and angles.

The extraction notes identify Section 26 of `cs3c00732_si_001.pdf` as the source of the published DFT coordinates. The manifest records page starts and Eopt in hartree; each XYZ comment records its reported optimized energy.
