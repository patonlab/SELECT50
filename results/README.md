# Core SELECT50 benchmark data

This folder contains [100 starting XYZ files](reference/starting_structures/), [energy tables](energetics/), and [geometry tables](geometry/). The [optimized geometry sets](../optimized-geometries/) are stored separately: eight comparison methods plus the reference method, each with 100 structures.

Energy table columns encode `minor_minus_major` and specify kcal/mol; the identity of `spe_corr` needs confirmation. Geometry tables use ångströms and bond indices are one-based. The bond table contains 1,512 coordinate-verified rows across nine methods. The supplied figure identifies `ja1g_minor` as P1–C3 bond breaking and `ja1g_major` as rotation; the latter has no assigned forming/breaking bond. See the [geometry notes](geometry/README.md) for assignments and the `ac32_minor` duplicate correction. Version labels are supplied labels, not independently verified checkpoints.

Case metadata, raw reference calculations, and the perturbed-guess study are missing. Major/minor names are provisional. See [repository status](../repository_status.md) and [contribution instructions](../CONTRIBUTING.md).
