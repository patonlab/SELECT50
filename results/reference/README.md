# Reference calculations and starting structures

[starting_structures/](starting_structures/) contains 100 XYZ files for the 50 SELECT50 pairs. Each file is identical to its optimized ωB97X-D3(BJ)/def2-TZVP reference copy. The electronic-energy reference is ωB97M-V/def2-TZVPD on those geometries.

## Files

| File | Contents |
| --- | --- |
| [case_metadata.csv](case_metadata.csv) | Case descriptions, formulas, selectivity classes, charges, multiplicities, metal classifications and source DOIs |
| [energy_orientation.csv](energy_orientation.csv) | 600 source/canonical orientation records: 50 cases × two energy quantities × six prediction methods |
| [failed_calculations.csv](failed_calculations.csv) | Four an82 geometry exclusions: both structures for g-xTB and r²SCAN-3c |
| [geometry_failed_calculations.csv](geometry_failed_calculations.csv) | Ten excluded measurement records associated with those four failed geometries |

## Roles and exclusions

Major is fixed by the lower reference electronic energy. ΔΔE and ΔΔG are calculated as minor minus major with the same structures across all methods. Five reference pairs reverse energetic preference between E and G. The orientation table records the conversion from source signs to these roles; 18 source pairs have an orientation factor of −1.

The an82 source labels map as `an82_decarb` → major and `an82_hat` → minor. Its g-xTB predictions are blank in both master energy tables, giving 49 evaluated pairs. Both ja1g structures are included in geometry/RMSD comparisons and excluded from bond analysis.

## Method labels

Core UMA models are UMA-M-1.1 and UMA-S-1.2. SPEcorr denotes high-level DFT single points on UMA-S geometries. Case metadata are repeated in the master energy tables; fields ending in `_inferred` are heuristic classifications.
