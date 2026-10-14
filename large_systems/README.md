# Large-system transition structures — Figure 7

This folder contains six optimized transition-state pairs (12 XYZ files), twelve matching DFT starting structures in [starting_DFT_structures/](starting_DFT_structures/), and [large_system_transition_pairs.csv](large_system_transition_pairs.csv). The study maps to manuscript Figure 7 and SI Table S22.

## Files, units and methods

The transition-pair CSV records the case, pathway, total atom count, charge, multiplicity, model label, electronic and Gibbs energies, imaginary frequency, selectivity class and recorded calculation time. Electronic/Gibbs energies are in hartree, frequencies in cm⁻¹, and `time_seconds` in seconds. Its model label is `uma_s_1p2` (UMA-S-1.2). XYZ coordinates use Å.

The manuscript figure compares each displayed major structure with its corresponding DFT starting geometry using an all-atom least-squares RMSD. Optimized and reference files are paired by filename and XYZ atom order.

## Case and file map

| Case | Optimized/reference filename stems | Literature DOI | Displayed major RMSD / Å |
| --- | --- | --- | ---: |
| Trichloromethylation | `trichloromethylation_major`, `trichloromethylation_minor` | `10.1039/C7OB01826F` | 0.273 |
| Hydroacetalization | `hydroacetalization_TS2`, `hydroacetalization_TS2prime` | `10.1021/jacs.4c04287` | 0.360 |
| Ir borylation | `borylation_C5`, `borylation_C6` | `10.1021/jacs.4c08027` | 0.331 |
| Prins-semipinacol | `Prins_semipinacol_Si`, `Prins_semipinacol_Re` | `10.1021/acscatal.4c01494` | 0.306 |
| Sulfide oxidation | `sulfide_oxidation_major`, `sulfide_oxidation_minor` | `10.1021/acs.orglett.4c02452` | 0.219 |
| β-Fluoride elimination | `beta_fluoride_elimination_TS3_major`, `beta_fluoride_elimination_TS4_minor` | `10.1021/jacs.1c05769` | 0.158 |

The six displayed-major RMSDs span 0.158–0.360 Å. Across all twelve structures, RMSDs span 0.158–1.090 Å under the same alignment convention.

In the case order above, the tabulated Gibbs differences round to 3.7, 3.1, 1.5, 1.7, 1.1 and 5.0 kcal/mol. The CSV's recorded major-job times span 234–714 s (3 min 54 s–11 min 54 s); all twelve recorded times span 213–826 s (3 min 33 s–13 min 46 s).

The figure's experimental comparisons are hydroacetalization 2.3, borylation >2.7, Prins-semipinacol 1.3 and sulfide oxidation 1.6 kcal/mol. The borylation value is a lower bound.
