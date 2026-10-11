# Geometry comparisons

[heavy_atom_rmsd.csv](heavy_atom_rmsd.csv) contains 800 heavy-atom RMSDs in ångströms. [ts_bond_distance.csv](ts_bond_distance.csv) contains 1,512 distances in ångströms: 168 assigned bonds across 99 structures, each evaluated for all nine geometry methods. Atom indices refer to the one-based order in each deposited XYZ file.

## Bond-table corrections

### ac32_minor

The r2SCAN-3c C17–H25 entry is **2.354753 Å**, matching [the deposited geometry](../../optimized-geometries/r2SCAN-3c/ac32_minor/ac32_minor.xyz). The conflicting duplicate reported 2.342934 Å and was removed. That value matches the separately retained UMA-S entry; the origin of its duplicate r2SCAN label is not established. All other existing rows were retained unchanged.

### ja1g

Assignments follow the supplied [render sheet 2](../../assets/figure_02.svg):

| Structure | Assignment | Bond-table treatment |
| --- | --- | --- |
| `ja1g_minor` | P1–C3 bond breaking | Added the P1–C3 distance calculated from each of the nine deposited method geometries. The reference distance is 2.383415 Å, consistent with the figure's 2.383 Å annotation. |
| `ja1g_major` | Rotation | No forming/breaking bond is assigned in the figure. This structure is intentionally excluded from the bond-distance table and explicitly accounted for in validation. |

Added distances were calculated as the Euclidean separation of atoms 1 and 3 and rounded to six decimal places, following the existing table. No coordinates were changed. The rotational assignment does not define a torsion's atom indices; no torsional coordinate has been inferred.
