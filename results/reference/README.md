# Reference calculations and starting structures

The [starting_structures/](starting_structures/) folder contains exactly 100 XYZ files, one `major` and one `minor` structure for each of 50 cases. They are reported as ωB97X-D3(BJ)/def2-TZVP geometries and are byte-identical to corresponding reference method files under `optimized-geometries/`.

## Draft case metadata and aliases

- [case_metadata.csv](case_metadata.csv) contains 50 rows transcribed from the working SI Table S2, keyed by the repository's stable `benchmark_id`. It includes formula, selectivity class, charge, multiplicity, source DOI, and the requested `mechanism_type` column. Mechanisms are proposals, not approved assignments.
- [structure_aliases.csv](structure_aliases.csv) contains 100 rows linking existing XYZ IDs to SI naming candidates and the original output labels recorded in the reference XYZ comments. It preserves the `ja37` ↔ `ja37f` case alias and legacy source names without renaming files.
- [metadata_notes.md](metadata_notes.md) defines fields, source provenance, validation results, and the remaining author decisions.

Major/minor roles remain provisional. Repository labels use reference electronic-energy ordering, whereas SI methods say Gibbs-energy ordering; five reference pairs reverse order between E and G. Every proposed manuscript role alias is marked `pending_U01`. Output filenames are evidence of local source labels, not verified original literature accessions.

Reference frequency and high-level single-point outputs are not deposited. The separately supplied SI documents a reference protocol, but original job settings/checkpoints are still needed to verify the deposited calculations. See [repository status](../../repository_status.md).
