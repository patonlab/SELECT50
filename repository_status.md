# Repository status

Updated 2026-10-12 against local HEAD `bc98c02021b8993215afd8d8235bdbc2b0a0c068` plus the working-tree integration described here. This is a local checkout audit; live GitHub verification was unavailable because credentials were missing. No commit, push or release was made.

## Current core inventory

- 50 cases including an82; 100 current reference copies, each byte-identical to the optimized DFT reference.
- 896 active optimized XYZs: seven methods × 100 and g-xTB/r²SCAN-3c × 98. Four obsolete ja30 copies are preserved in the historical archive.
- Two 50-row energy masters; g-xTB has 49 evaluated pairs in each. an82 source description/DOI are filled from current SI S2.
- 796 heavy-atom RMSDs, 1,488 coordinate-verified bond distances and 1,322 bond errors.
- 166 reference bonds over 98 structures. Both ja1g structures are excluded by author confirmation on 2026-10-12. g-xTB/r²SCAN-3c each cover 163 bonds over 96 structures.
- 600 source-orientation records, case metadata and explicit geometry-failure ledgers. Supporting aliases, method/provenance notes and source hashes remain in the Downloads audit/export at the author’s request.

Major is fixed by lower reference electronic energy and retained for Gibbs energy. The five negative reference G gaps are valid. Eighteen original pairs were jointly reoriented. Model labels for the core are UMA-M-1.1 and UMA-S-1.2; these do not establish per-job checkpoint identity.

## Validation

Run `python3 scripts/validate_data.py` and `python3 -m unittest discover -s tests` (install NumPy to run the RMSD tests). The validator covers XYZ parsing, current method/case coverage, reference-copy equality, element ordering, metadata atom counts, required CSV columns, allowed energy missingness, common E/G orientations, energy method summaries, bond distances/errors and local inline Markdown links. Regression tests reject coverage gaps, undeclared missing values and independent G-role reversal.

The RMSD calculator separately reproduces all 796 available values at six-decimal output precision. Its strict mode reports the four known absent geometries as failures. Passing structural validation does not verify every original job, plot, category statistic, literature assignment or manuscript statement.

## Other deposited work and remaining gaps

Recovery includes 300 logs and 299 final XYZ exports plus a placeholder for the documented missing jo8v final export. Reproducible recovery summaries and generator/reference provenance remain pending. Large systems contain 12 optimized and 12 starting geometries plus E/G/frequency/timing values; original job versions, mappings and timing scope need reconciliation.

Cycloaddition and palladium landscape data are deposited. C05/C06 corrections, source-backed audits and replacement plots are prepared in the user's separate Downloads audit and **have not been integrated here**, as requested. Existing study guides and manuscript/SI still need those corrections.

Calibration supplies 13 Ni and 15 Pd ee entries and saved 100-partition analyses. Pd experimental P03/P11 values conflict across current inputs, and the Pd script requests 1000 partitions. Resolve inputs/protocol before regenerating final results. Calibration OLE drawings still need portable replacements.

The current manuscript/SI remain unchanged: core tables/figures and study corrections still require a careful merge. Source identity questions are tracked in [reference guide](results/reference/README.md). Original reference outputs, complete checkpoint/environment/thermal provenance, final release attribution and an actual archival DOI remain outstanding. Input templates are examples, not original-job evidence.

## Historical audit results

The pre-integration validator failed because it expected obsolete folders, CSV schemas and complete coverage. Eight broken links and the stale reference/MACE ja30 copies were found then. Those findings describe the earlier working tree; use the current validator and inventory above for this integration. C05/C06 scientific findings remain staged for review.
