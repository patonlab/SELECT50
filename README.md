# SELECT50

Data and benchmark materials for *Mapping Organic Selectivity with Machine-Learned Interatomic Potentials: Introduction of the SELECT50 Data Set*.

SELECT50 contains 50 selectivity cases and 100 reference transition-state geometries, reported as optimized at ωB97X-D3(BJ)/def2-TZVP. Structure names use `major` and `minor` according to the reference electronic-energy ordering. These labels remain provisional until the author's final naming correction is supplied.

## Available data

- [100 starting XYZ files](results/reference/starting_structures/), named `<case_id>_major.xyz` and `<case_id>_minor.xyz`.
- [Optimized geometries](optimized-geometries/): eight comparison methods and the reference method, each with 100 XYZ files. Reference method files are byte-identical to the starting structures.
- [ΔE](results/energetics/delta_e.csv) and [ΔG](results/energetics/delta_g.csv): 50 cases each, with differences expressed as `minor_minus_major`. Energy units and calculation settings are not documented in the supplied files; confirm them before interpreting values.
- [Heavy-atom RMSD](results/geometry/heavy_atom_rmsd.csv): 800 rows, in ångströms.
- [TS bond distances](results/geometry/ts_bond_distance.csv): 1,504 rows, in ångströms, with one-based atom indices. Main rows cover 98 structures; the `ja1g` pair is missing. One auxiliary row has an unresolved method label; see [repository status](repository_status.md).
- [Reaction landscapes](reaction_landscapes/): multi-frame XYZ files containing 10 cycloaddition structures and 21 palladium-cyclization structures.

## Repository layout

```text
optimized-geometries/<method>/<structure_id>/<structure_id>.xyz
results/
  reference/starting_structures/    100 benchmark input XYZ files
  energetics/                      ΔE and ΔG tables
  geometry/                        RMSD and TS bond-distance tables
  perturbed_guesses/               placeholder
reaction_landscapes/               two multi-frame XYZ files
large_systems/                     placeholder
calibration/                       two screening-study placeholders
scripts/                          data validation
```

The manuscript, Supporting Information, naming maps, raw calculation outputs, and model checkpoints are not included. The [status file](repository_status.md) records outstanding data and metadata needs.

## License and citation

Scientific data and documentation are licensed under [CC BY 4.0](LICENSE-DATA). Source code, tests, and repository automation and configuration are licensed under [MIT](LICENSE-CODE). See [LICENSE](LICENSE) for the scope of each license.

Use [CITATION.cff](CITATION.cff) or GitHub's **Cite this repository** feature to cite the repository. Include the commit hash used to identify the exact snapshot. A plain-text citation is:

> Pounder, Austin. SELECT50: Data and benchmark materials. https://github.com/apounder/SELECT50. Commit: <commit hash used>.

This citation covers the repository. The associated manuscript's complete author list and publication DOI have not been supplied.

## Validate the data

Run with Python 3.10 or newer; no third-party packages are required:

```sh
python3 scripts/validate_data.py
python3 -m unittest discover -s tests
```

The validator checks XYZ frames, case and method coverage, atom ordering, CSV schemas, finite numeric values, heavy-atom counts, bond indices and symbols, bond distances against coordinates, and local Markdown links. It reports the known auxiliary bond row and missing `ja1g` bond definitions separately. It does not establish calculation provenance or reproduce energy or RMSD calculations.

## Contribute a benchmark

See [CONTRIBUTING.md](CONTRIBUTING.md) and the [data guide](results/README.md). Submissions should identify model and version, input coordinates, settings, raw outputs, optimized XYZ files, energy units, and failed cases. Retain current case IDs and pair ordering until the naming map is finalized.
