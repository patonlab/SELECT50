# SELECT50

![SELECT50 banner showing two selectivity energy profiles](assets/select50-github-banner.png)

**Molecular geometries and energetics for studying organic and organometallic reaction selectivity with machine-learned interatomic potentials.**

SELECT50 examines competing transition-state structures, their relative electronic and Gibbs energies, and the structural changes obtained with different computational methods. The repository accompanies *Mapping Organic Selectivity with Machine-Learned Interatomic Potentials: Introduction of the SELECT50 Data Set*.

The data cover transition-state geometry comparisons, selectivity-energy benchmarks, recovery from perturbed starting geometries, reaction landscapes, larger molecular systems, and calibration against experimental selectivity.

## Transition-state gallery

The four sheets illustrate 50 transition-state pairs. Each pair is labelled by case ID, with major and minor structures shown side by side. Select a sheet to open the full-size PNG.

[![SELECT50 transition-state gallery, sheet 1](assets/figure_01.png)](assets/figure_01.png)

[![SELECT50 transition-state gallery, sheet 2](assets/figure_02.png)](assets/figure_02.png)

[![SELECT50 transition-state gallery, sheet 3](assets/figure_03.png)](assets/figure_03.png)

[![SELECT50 transition-state gallery, sheet 4](assets/figure_04.png)](assets/figure_04.png)

## Repository contents

| Dataset | Contents | Location |
| --- | --- | --- |
| Reference and optimized geometries | Reference coordinates and method-specific transition-state XYZ files | [Reference structures](results/reference/starting_structures/), [optimized geometries](optimized-geometries/) |
| Selectivity energies | Electronic- and Gibbs-energy separation tables, each containing 50 case rows | [Electronic energies](results/energetics/delta_delta_e.csv), [Gibbs energies](results/energetics/delta_delta_g.csv) |
| Geometry comparisons | 796 heavy-atom RMSD rows and 1,488 bond-distance rows | [RMSDs](results/geometry/heavy_atom_rmsd.csv), [bond distances](results/geometry/ts_bond_distances.csv) |
| Transition-state recovery | 300 calculation logs and 300 XYZ exports, organized by transition structure, perturbation magnitude, and replicate | [Recovery calculations](ts_recovery/), [summary table](results/perturbed_guesses/ts_recovery_full_data.csv) |
| Cycloaddition landscape | Ten stationary structures with E/G records, a 100-point two-dimensional scan, IRC trajectories, and 72 activation-strain rows | [Cycloaddition data](reaction_landscapes/cycloaddition/) |
| Palladium-cyclization landscape | 20 UMA stationary structures, spectator-adjusted energy profiles, IRC and NEB paths, a 73-point dihedral scan, and 21 published DFT structures | [Palladium data](reaction_landscapes/palladium_cyclization/) |
| Larger systems | Six transition-state pairs: 12 optimized structures, 12 DFT starting structures, and a table of energies, frequencies, and recorded timings | [Structures](large_systems/), [transition-pair table](large_systems/large_system_transition_pairs.csv) |
| Experimental calibration | Ni and Pd series with 13 and 15 entries, respectively, plus transition structures, fitting scripts, saved split memberships, and held-out predictions | [Nickel series](calibration/nickel_epoxide_arylation/), [palladium series](calibration/palladium_enone_addition/) |

## Computational methods

The optimized-geometry directories represent the following methods:

| Category | Methods |
| --- | --- |
| Reference DFT | ωB97X-D3(BJ)/def2-TZVP |
| Machine-learned interatomic potentials | eSEN-medium, MACE-OMOL, UMA-M-1.1, UMA-S-1.2 |
| Additional DFT methods | PBE0-D3(BJ)/def2-SVP, r²SCAN-3c, ωB97X-D3(BJ)/def2-SVP |
| Semiempirical method | g-xTB |

The recovery logs identify UMA-S-1.2 with the OMol task. Reaction-landscape records identify UMA-S-1.1 with the OMol task through ORCA's external-method interface.

The active core contains 50 pairs, including an82. Seven geometry sets contain 100 structures each; g-xTB and r²SCAN-3c contain 98 each because both an82 calculations failed. The g-xTB energy tables evaluate 49 pairs. Both ja1g structures are excluded from bond analysis, as confirmed by the author on 2026-10-12. See [metadata and exclusions](results/reference/README.md) and [current repository status](repository_status.md).

## Data conventions

- **Coordinates and geometry:** XYZ coordinates, RMSDs, and bond distances use ångströms (Å). Bond-table atom indices are one-based.
- **Core energy differences:** CSV headers specify kcal/mol. ΔΔE = E(minor) − E(major), with the analogous definition for ΔΔG. Major is the lower reference electronic-energy structure; that role is fixed for both E and G. Only reference ΔΔE is positive by construction; five reference ΔΔG values are negative.
- **Absolute energies:** Landscape and large-system tables report electronic and Gibbs energies in hartree. Relative landscape energies are recorded in kcal/mol.
- **Frequencies and timings:** Imaginary frequencies in the large-system table use cm⁻¹; its `time_seconds` field uses seconds.
- **Structure names:** Core XYZ filenames use `<case_id>_major.xyz` and `<case_id>_minor.xyz`. Method geometries are organized as `optimized-geometries/<method>/<structure_id>/<structure_id>.xyz`.
- **Heavy-atom RMSD:** The supplied calculation script pairs atoms by XYZ order, excludes hydrogen, and uses an equally weighted Kabsch alignment with proper rotation and no outlier rejection.

Recovery runs are arranged at six target Cartesian RMSD magnitudes—0.05, 0.08, 0.11, 0.14, 0.17, and 0.20 Å—with five replicates for each of ten transition structures. Their logs include optimization status, step counts, electronic energies, vibrational frequencies, and thermochemical corrections.

The saved Ni and Pd calibration analyses each contain 100 partitions. Their tables record training/test memberships, fitted coefficients, held-out predictions, and raw and calibrated error metrics.

## Access and analysis

Clone the repository and record the commit used for an analysis:

```sh
git clone https://github.com/apounder/SELECT50.git
cd SELECT50
git rev-parse HEAD
```

Coordinates are distributed as XYZ files, tabulated results as CSV files, and recovery output as text logs. These formats can be read directly with molecular visualization tools, spreadsheets, or analysis scripts.

The repository includes:

- A [MACE-OMOL/ASE calculation example](examples/calculations/opt_mace_omol.py) with setup and run instructions.
- ORCA external-UMA templates for [minimum optimization](examples/calculations/job_orca.inp), [TS optimization](examples/calculations/job_orca_ts.inp), [relaxed scans](examples/calculations/job_orca_scan.inp), and [IRC calculations](examples/calculations/job_orca_irc.inp).
- [A structural data validator](scripts/validate_data.py) for XYZ parsing, table schemas, case coverage, atom ordering, bond distances, and local documentation links.
- [A heavy-atom RMSD calculator](scripts/calculate_heavy_atom_rmsd.py) using NumPy and the ωB97X-D3(BJ)/def2-TZVP geometry directory as its reference.
- [Ni calibration](calibration/nickel_epoxide_arylation/calibrate.py) and [Pd calibration](calibration/palladium_enone_addition/calibrate.py) scripts using NumPy for repeated-split linear fits of experimental ΔΔG against computed ΔΔG.
- [Python tests](tests/) for the geometry-analysis and validation code.

## Citation

Use the repository's [citation metadata](CITATION.cff) when citing these data or scripts. Include the repository URL and the commit hash used in your analysis to identify the data version.

## License

Scientific data and documentation are licensed under [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE-DATA). Source code, tests, and repository automation and configuration are licensed under the [MIT License](LICENSE-CODE). See [LICENSE](LICENSE) for the scope of each license. Separately identified third-party material retains its original license and attribution notices.

## Contributions

Corrections and contributions are welcome through GitHub issues and pull requests. See [CONTRIBUTING.md](CONTRIBUTING.md) for contribution guidelines.
