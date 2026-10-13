# Large-system transition structures — Figure 8

This folder contains all six optimized TS pairs (12 XYZ files) and twelve matching filenames under [starting_DFT_structures/](starting_DFT_structures/). XYZ coordinates and RMSDs are in Å. The current figure reports all-atom least-squares RMSDs for the displayed major-product structures; it does not report a pair-averaged heavy-atom RMSD.

The separately supplied working manuscript/SI map these data to Figure 8 and Table S22. Those DOCX drafts are available beside the checkout for review, not deposited in this repository. S22 contains E/G in hartree and imaginary frequencies in cm⁻¹ for all twelve TSs; a deposited energy/frequency table has not yet been exported.

## Case and reference map

References and RMSDs below are transcribed from the supplied drafts, not recalculated or independently checked against the original publications.

| Case | Optimized/reference filename stems | Main reference | SI reference | DOI | Displayed RMSD / Å |
| --- | --- | ---: | ---: | --- | ---: |
| Trichloromethylation | `trichloromethylation_major`, `trichloromethylation_minor` | 85 | 81 | `10.1039/C7OB01826F` | 0.273 |
| Hydroacetalization | `hydroacetalization_TS2`, `hydroacetalization_TS2prime` | 86 | 82 | `10.1021/jacs.4c04287` | 0.360 |
| Ir borylation | `borylation_C5`, `borylation_C6` | 12 | 83 | `10.1021/jacs.4c08027` | 0.331 |
| Prins-semipinacol | `Prins_semipinacol_Si`, `Prins_semipinacol_Re` | 84 | 84 | `10.1021/acscatal.4c01494` | 0.655 |
| Sulfide oxidation | `sulfide_oxidation_major`, `sulfide_oxidation_minor` | 87 | 85 | `10.1021/acs.orglett.4c02452` | 0.219 |
| β-Fluoride elimination | `beta_fluoride_elimination_TS3_major`, `beta_fluoride_elimination_TS4_minor` | 88 | 86 | `10.1021/jacs.1c05769` | 0.158 |

The two XYZ sets have matching atom counts and element sequences for each filename. This does not independently prove correct correspondence of identical-element atoms or authenticate the original published structures. Original structure labels/accessions and ONIOM/full-system extraction details still need a source map.

S22's paired Gibbs differences round to 3.7, 3.1, 1.5, 1.7, 1.1, and 5.0 kcal/mol in the row order above. Explicit filename-to-S22 major/minor assignments are still needed for TS2/TS2prime and C5/C6. Prins Si/Re electronic energies in the XYZ comments correspond to S22's major/minor rows, respectively.

## Methods and remaining records

Model versions need reconciliation. Figure 8/S22/large-system prose say UMA-S-1.1; both Prins XYZ comments explicitly identify UMA Small 1.2, fairchem 2, and Rowan. The top-level README previously generalized all non-profile UMA-S jobs as 1.2. Original checkpoint/job records must establish the intended per-case versions before those labels or data are changed.

The SI describes ASE/geomeTRIC TS optimization and harmonic free energies at 300 K for this study. It does not replace the missing per-job calculation inputs, frequency exports, atom-map/alignment settings, or complete provenance.

Experimental values are only partially reported in Figure 8: hydroacetalization 2.3, borylation >2.7, Prins 1.3, and sulfide oxidation 1.6 kcal/mol. No parenthetical experimental number is given for trichloromethylation or β-fluoride elimination. Deposit raw ee/er or regioselectivity ratio, actual temperature, source location, sign convention, and exact/bound status per case. The borylation lower bound should remain a bound in error comparisons.

The main text reports RTX 5070 and a timing range of 3 min 54 s–11 min 54 s, but a per-case timing table/log map is not deposited. Define whether each time refers to one TS or a pair and which optimization/frequency tasks it includes. Exporting S22, adding its RMSD column, reproducing overlays, and completing experimental/timing metadata remain later Figure 8 tasks.
