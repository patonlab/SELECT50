![SELECT50 banner showing two selectivity energy profiles](select50-github-banner.png)

# SELECT50

**A benchmark dataset for evaluating machine-learned interatomic potentials in organic selectivity.**

This repository contains molecular geometries, energy differences, geometry comparisons, and supporting case studies for *Mapping Organic Selectivity with Machine-Learned Interatomic Potentials: Introduction of the SELECT50 Data Set*.

The core benchmark comprises **50 selectivity cases and 100 reference transition-state geometries**, reported as optimized at ωB97X-D3(BJ)/def2-TZVP. Eight comparison methods are represented alongside the reference geometries. Supporting datasets cover reaction landscapes, larger molecular systems, and calibration against experimental selectivity.

**Release status:** This is a research snapshot. Some experimental metadata, computational provenance, and figure-reproduction materials remain incomplete. A publication DOI and archived dataset DOI are not currently recorded in this repository.

## Data inventory

| Dataset | Contents | Location |
| --- | --- | --- |
| Reference structures | 100 XYZ files: two transition states per selectivity case | [Starting structures](results/reference/starting_structures/) |
| Optimized structures | Nine geometry sets, each containing 100 XYZ files; reference copies are byte-identical to the starting structures | [Optimized geometries](optimized-geometries/) |
| Selectivity energies | Electronic-energy and Gibbs-energy differences for 50 cases | [ΔE](results/energetics/delta_e.csv), [ΔG](results/energetics/delta_g.csv) |
| Geometry comparisons | 800 heavy-atom RMSDs and 1,504 bond-distance rows | [Geometry tables](results/geometry/) |
| Cycloaddition landscape · Figure 6 | Ten stationary-point structures, Gibbs energies, a 100-point two-dimensional scan, IRC trajectories, and activation-strain tables | [Cycloaddition](reaction_landscapes/cycloaddition/) |
| Palladium-cyclization landscape · Figure 7 | 21 UMA stationary-point structures, Gibbs energies, spectator-adjusted profiles, IRC/NEB paths, a 73-point dihedral scan, and published DFT geometries | [Palladium cyclization](reaction_landscapes/palladium_cyclization/) |
| Larger systems · Figure 8 | Six transition-state pairs, comprising 12 XYZ files | [Large systems](large_systems/) |
| Selectivity calibration · Figure 9 | Ni and Pd transition-state pairs, calibration scripts, split memberships, held-out predictions, and fit/error summaries | [Ni study](calibration/nickel_epoxide_arylation/), [Pd study](calibration/palladium_enone_addition/) |
| Perturbed-guess recovery · Figure 5 | Reserved directory; recovery inputs and results are not yet deposited | [Perturbed guesses](results/perturbed_guesses/) |

## Data conventions

- **Coordinates and distances:** XYZ coordinates, RMSDs, and bond distances are expressed in ångströms (Å). Bond-table atom indices are one-based.
- **Energy differences:** The core ΔE and ΔG headers specify kcal/mol. Differences use `minor_minus_major`. Reaction-landscape tables distinguish absolute energies in hartree from relative energies in kcal/mol where indicated in their headers.
- **Structure identifiers:** Core XYZ files use `<case_id>_major.xyz` and `<case_id>_minor.xyz`. Major/minor labels follow the supplied reference electronic-energy ordering and remain provisional pending the final naming map.
- **Methods:** Geometry folders include eSEN, g-xTB, MACE-OMOL, UMA-M, UMA-S, PBE0-D3(BJ)/def2-SVP, r2SCAN-3c, and ωB97X-D3(BJ)/def2-SVP, alongside the ωB97X-D3(BJ)/def2-TZVP reference. Folder names identify the supplied comparison sets; checkpoint and software provenance require separate verification.
- **Thermochemistry and experimental selectivity:** Use the conditions and definitions associated with each study. The meaning of `spe_corr`, the full core thermochemistry protocol, and experimental selectivity metadata are not yet documented completely.

## Repository structure

```text
SELECT50/
├── optimized-geometries/<method>/<structure_id>/<structure_id>.xyz
├── results/
│   ├── reference/starting_structures/
│   ├── energetics/
│   ├── geometry/
│   └── perturbed_guesses/
├── reaction_landscapes/
│   ├── cycloaddition/
│   └── palladium_cyclization/
├── large_systems/
├── calibration/
│   ├── nickel_epoxide_arylation/
│   └── palladium_enone_addition/
├── scripts/
├── tests/
├── CITATION.cff
└── LICENSE, LICENSE-DATA, LICENSE-CODE
```

## Access and use

Clone the repository or download an archive from GitHub:

```sh
git clone https://github.com/apounder/SELECT50.git
cd SELECT50
git rev-parse HEAD
```

Record the commit hash with any analysis to identify the exact dataset snapshot. XYZ and CSV files can be read without installing project-specific software. Consult the [data guide](results/README.md) for the core benchmark and retain the supplied case IDs when comparing methods.

### Validation

The validation tools require Python 3.10 or newer and use only the standard library:

```sh
python3 scripts/validate_data.py
python3 -m unittest discover -s tests
```

The validator checks XYZ syntax, core case/method coverage, atom ordering, CSV schemas, bond distances, and local Markdown links. It currently expects the previous palladium-trajectory location and energy-column names, so it does **not** pass on this snapshot. The seven existing regression tests pass, but do not establish complete dataset validity or reproduce the reported scientific results.

### Calibration analyses

The calibration scripts require NumPy and implement repeated-split linear fits of experimental ΔΔG against raw predicted ΔΔG. Saved outputs include train/test memberships, fitted coefficients, and held-out errors.

Independent reruns currently require input repairs: the Ni series headers differ from those expected by its script, and `pd_series.csv` is not deposited. See the [Ni analysis notes](calibration/nickel_epoxide_arylation/README.txt) for the seed, split protocol, and additional provenance limitations.

## Reproducibility and outstanding metadata

The deposited files support inspection of geometries and numerical comparisons. The following items remain necessary for a complete publication package:

- Figure 5 perturbation inputs or their generator, recovery results, bootstrap procedure, and consistent acceptance criteria.
- Per-case experimental selectivity, temperature, literature references, large-system reference geometries/RMSDs, and timings.
- Complete calibration input tables and identification of the analyses used in the manuscript and SI figures.
- Computational settings, model/checkpoint versions, source outputs, and published comparison data sufficient to reproduce the reported calculations and overlays.
- Resolution of the duplicate `ac32_minor` r2SCAN-3c bond entry and the missing `ja1g` bond definitions; clarification of palladium NEB image counts and geometry-comparison conventions.
- A case-metadata table with mechanism types, finalized naming, updated validation/documentation, and a versioned archival deposit.

The manuscript, its Supporting Information, and model checkpoints are not deposited here. The [repository status](repository_status.md) includes an earlier audit; its inventory and some unresolved-item descriptions predate subsequent data additions.

## Citation

Please cite this repository when using its data or analysis scripts, and identify the commit or archived version used. Machine-readable citation metadata is provided in [CITATION.cff](CITATION.cff), which also enables GitHub's **Cite this repository** feature.

> Pounder, Austin. *SELECT50: Data and benchmark materials*. https://github.com/apounder/SELECT50.

The current citation identifies the repository. It does not substitute for the associated manuscript's complete author list and publication reference, which will need to be recorded when available.

## License

Scientific data and documentation are licensed under [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE-DATA). Source code, tests, and repository automation and configuration are licensed under the [MIT License](LICENSE-CODE). See [LICENSE](LICENSE) for the scope of each license. Separately identified third-party material retains its original license and attribution notices.

## Contributions and questions

Corrections, additional calculation provenance, and benchmark submissions are welcome through GitHub issues and pull requests. See [CONTRIBUTING.md](CONTRIBUTING.md) for submission requirements. Include method/version identifiers, input coordinates, computational settings, energy units, optimized structures, raw outputs, and failed cases. Retain current case IDs and pair ordering until the naming map is finalized.
