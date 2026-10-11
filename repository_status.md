# SELECT50 repository status

Audited 2026-10-09 against files available in this workspace. This data repository is intended for local analysis and benchmark submissions; there is no application server.

## Original audit inventory (2026-10-09)

| Data | Available |
| --- | --- |
| Reference inputs | 100 XYZ files, 50 complete major/minor pairs |
| Optimized geometries | Nine method folders, 100 XYZ files each; eight comparison methods and the reference method |
| Reference copies | All 100 reference-method files are byte-identical to benchmark inputs |
| Electronic-energy differences | 50 unique cases |
| Free-energy differences | 50 unique cases |
| Heavy-atom RMSD | 800 rows: 100 structures and eight methods |
| TS bond distances | 1,503 main rows: 98 structures and nine methods, plus one auxiliary row |
| Cycloaddition landscape | 10 XYZ frames |
| Palladium-cyclization landscape | 21 XYZ frames |
| Perturbed guesses, large systems, calibration | Empty placeholders |

All 1,002 XYZ files parse (1,031 frames total). Every optimized set matches reference IDs and atom ordering. RMSD heavy-atom counts match reference inputs. All 1,503 main bond-distance rows match method coordinates within 0.00002 Å. Reference ΔE values are positive under the supplied `minor_minus_major` convention. Scientific file bytes were preserved during this audit.

## Repository fixes

- Corrected `reults/` to [results/](results/).
- Replaced documentation for absent `select50/`, `paper/`, and `provenance/` trees with actual links and counts.
- Removed unsupported claims that naming maps, a move manifest, manuscripts, and 30 unassigned geometries were included and verified.
- Added a standard-library validator and regression tests. Automatic GitHub Actions validation was subsequently removed at the repository owner's request; checks can still be run manually.
- Added ignore rules for local credentials, caches, editor files, and temporary files.

## Unresolved scientific metadata

1. **Energy units and provenance:** Neither energy CSV declares units. Source calculations, thermochemistry settings, and the meaning of `spe_corr` are missing. Confirm these before interpreting numerical comparisons; do not guess units or silently convert values.
2. **Bond-table corrections (2026-10-10):** Removed the conflicting `ac32_minor` r2SCAN-3c C17–H25 duplicate, retaining the coordinate-matched 2.354753 Å row. Added the figure-assigned `ja1g_minor` P1–C3 breaking bond for all nine methods. The figure identifies `ja1g_major` as rotation with no assigned forming/breaking bond, so it is intentionally excluded. The corrected table contains 1,512 verified distances covering 99 structures. See [geometry notes](results/geometry/README.md) for provenance. The rotational torsion indices remain unspecified.
3. **Naming and versions:** The corrected major/minor map is absent. Earlier notes flag `ja37` versus `ja37f`, legacy manuscript signs (including `jo2y`), and UMA 1.1 versus 1.2 labels. The manuscript and checkpoints are unavailable, so these notes remain unverified.
4. **Reproducibility:** Raw outputs, settings, energies, and checkpoints are incomplete. Bond distances can be checked against coordinates. The supplied [RMSD script](scripts/calculate_heavy_atom_rmsd.py) uses NumPy Kabsch fitting, excludes hydrogen, and pairs atoms by XYZ order without outlier rejection. All 800 core RMSDs were reproduced within 0.000001 Å using Python 3.14.4 and NumPy 2.5.4. The [geometry methods](results/geometry/README.md#rmsd-methodology) document the calculation and distinguish it from the separately reported PyMOL `align` protocol. The original analysis environment and provenance of other supporting RMSD tables remain incomplete.
5. **Other studies:** Landscape coordinates are present, but raw outputs and path-analysis settings are missing. Perturbed guesses, large-system pairs, and calibration data remain absent. Earlier notes mention six large-system pairs and 300 perturbed runs; neither count can be verified.
6. **Publication metadata:** Manuscript, Supporting Information, the manuscript's complete author list, and publication DOI are absent. Repository citation metadata is now available in [CITATION.cff](CITATION.cff). The repository owner selected CC BY 4.0 for data and documentation and MIT for source code and repository automation and configuration; see [LICENSE](LICENSE).

These limitations allow a clearly labeled preliminary snapshot to be committed, but prevent presenting it as a complete, reproducible publication package.

## Git readiness

Checked 2026-10-10: the workspace is an initialized Git repository on `main`, tracking `origin/main` at `https://github.com/apounder/SELECT50.git`, with two existing commits. Existing history and working-tree data changes were preserved. Local credentials and agent configuration are ignored. Text checkout rules use LF; scientific XYZ and CSV files retain their supplied bytes. No commit or push was performed for the licensing and citation update.

## Repeatable checks

```sh
python3 scripts/validate_data.py
python3 -m unittest discover -s tests
```

Validation establishes file consistency, not scientific correctness or calculation provenance.

The audit run passed the full validator and seven regression tests. SHA-256 checks confirmed that all 1,006 scientific files (1,002 XYZ and four CSV files) remain byte-identical after relocation. Text whitespace checks passed. A working-tree scan found no matches for common AWS access keys, GitHub tokens, OpenAI-style API keys, or private-key blocks; this scan does not cover unavailable Git history or every possible secret format. These results describe the original audit snapshot. The automatic GitHub Actions validation workflow has since been removed.
