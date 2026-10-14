# Reference calculations and starting structures

[starting_structures/](starting_structures/) contains 100 XYZ files for the current 50 cases, including an82. Every file is byte-identical to its ωB97X-D3(BJ)/def2-TZVP optimized reference copy.

- [case_metadata.csv](case_metadata.csv): 50 cases, formula, selectivity class, charge, multiplicity and source identifiers.
- [energy_orientation.csv](energy_orientation.csv): 600 source/canonical orientation records (50 cases × two quantities × six prediction methods). E and G use the same roles; 18 source pairs are reoriented.
- [failed_calculations.csv](failed_calculations.csv): four an82 geometry exclusions, covering both structures for g-xTB and r²SCAN-3c.
- [geometry_failed_calculations.csv](geometry_failed_calculations.csv): source ledger of ten removed measurement rows arising from those four failed geometries; these are not ten independent calculation failures.

Case metadata are repeated in the electronic- and Gibbs-energy tables. The ac16 source DOI is `10.1021/acscatal.1c00616`; the ac30 source DOI is `10.1021/acsomega.9b00030`.

## Roles and exclusions

Major is fixed by lower reference electronic energy. The same role applies to Gibbs energies; five reference pairs reverse energetic preference between E and G. an82 source names are `an82_decarb` → major and `an82_hat` → minor. Both g-xTB an82 energy predictions remain blank, giving 49 evaluated pairs. The validator enforces these two energy exclusions directly; they are not additional rows in the geometry-failure ledger.

Both ja1g structures are excluded from bond analysis and included in the geometry/RMSD datasets.

## Computational labels

The core's supplied checkpoint labels are UMA-M-1.1 and UMA-S-1.2. The electronic reference is ωB97M-V/def2-TZVPD on ωB97X-D3(BJ)/def2-TZVP geometries. SPEcorr denotes high-level DFT single points on UMA-S geometries.

Fields ending in `_inferred` in the wide energy tables are heuristic classifications. The [validator](../../scripts/validate_data.py) checks metadata agreement, pair coverage, retained E/G orientations and declared exclusions.
