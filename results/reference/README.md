# Reference calculations and starting structures

[starting_structures/](starting_structures/) contains 100 XYZ files for the current 50 cases, including an82. Every file is byte-identical to its ωB97X-D3(BJ)/def2-TZVP optimized reference copy. Four obsolete ja30 reference/MACE files were moved intact to the [historical archive](../../archive/replaced_ja30/README.md).

- [case_metadata.csv](case_metadata.csv): 50 cases, formula, selectivity class, charge, multiplicity and source identifiers. an82's description and DOI were completed from the current supplied SI Table S2; no numerical data changed.
- [energy_orientation.csv](energy_orientation.csv): 600 source/canonical orientation records (50 cases × two quantities × six prediction methods). E and G use the same roles; 18 source pairs are reoriented.
- [failed_calculations.csv](failed_calculations.csv): four an82 geometry exclusions, covering both structures for g-xTB and r²SCAN-3c.
- [geometry_failed_calculations.csv](geometry_failed_calculations.csv): source ledger of ten removed measurement rows arising from those four failed geometries; these are not ten independent calculation failures.

Metadata and orientation records were imported from the supplied SELECT50_csv_metadata export. Supporting import hashes, aliases, method/provenance notes and source-review records remain in the separate Downloads audit/export rather than being added here, at the author's request.

## Roles and exclusions

Major is fixed by lower reference electronic energy. The same role applies to Gibbs energies; five reference pairs reverse energetic preference between E and G. an82 source names are `an82_decarb` → major and `an82_hat` → minor. Both g-xTB an82 energy predictions remain blank, giving 49 evaluated pairs. The validator enforces these two energy exclusions directly; they are not additional rows in the geometry-failure ledger.

On 2026-10-12 the author confirmed that both ja1g structures remain excluded from bond analysis, matching the current SI and tables. Both remain in the geometry/RMSD datasets. This supersedes the old minor P1–C3 inclusion, without changing numerical tables.

## Provenance limits and remaining review

The core's supplied checkpoint labels are UMA-M-1.1 and UMA-S-1.2; names are not per-job file hashes. The reported electronic reference is ωB97M-V/def2-TZVPD on ωB97X-D3(BJ)/def2-TZVP geometries. SPEcorr uses high-level single points on UMA-S geometries; complete original thermal/job provenance remains outstanding.

Fields ending in `_inferred` in the wide energy tables are heuristic classifications, not approved mechanism assignments or a four-family census. Filling an82's description does not validate or regenerate those fields. Source review remains open for the ac30 incomplete DOI, ac16 DOI conflict, ja37/ja37f identity, competing-spin cases and ja1m/jo90/sc99 assignments. Matching element sequences across geometries does not prove correspondence among identical atoms.

Reference frequency and high-level single-point outputs are not deposited. Example input files describe related workflows but do not establish settings of every original job. See [repository status](../../repository_status.md).
