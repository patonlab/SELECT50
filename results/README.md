# Core SELECT50 benchmark data

This folder contains [100 starting XYZ files](reference/starting_structures/), [energy tables](energetics/), and [geometry tables](geometry/). The [optimized geometry sets](../optimized-geometries/) are stored separately: eight comparison methods plus the reference method, each with 100 structures.

Energy table columns encode `minor_minus_major`; their units and the identity of `spe_corr` need confirmation. Geometry tables use ångströms and bond indices are one-based. The bond table omits both `ja1g` structures and includes one unresolved auxiliary row. Version labels are supplied labels, not independently verified checkpoints.

Case metadata, raw reference calculations, and the perturbed-guess study are missing. Major/minor names are provisional. See [repository status](../repository_status.md) and [contribution instructions](../CONTRIBUTING.md).
