# Cycloaddition reaction landscape

This directory contains the [ten stationary-point geometries](cycloaddition_structures.xyz), [electronic and Gibbs energies](gibbs_energies.csv), [two-dimensional scan](2d-scan/scan2d.csv), [IRC trajectories and energy tables](irc/), and [activation-strain decomposition](activation_strain.csv).

This study maps to manuscript Figure 5, SI Tables S19/S20, and Figure S13. The literature source DOI is `10.1021/acs.orglett.6c01569`.

## Figure inputs and units

- Figure 5B/S19: ten stationary-point E/G rows in `gibbs_energies.csv`; absolute energies in hartree, relative energies in kcal/mol, recorded temperature 298.15 K, pressure 1 atm, quasi-RRHO flag, and imaginary frequencies in cm⁻¹.
- Figure 5C: 100 points in `2d-scan/scan2d.csv`, headers `r1,r2,E_rel_kcal`; distances and XYZ coordinates in Å, relative electronic energies in kcal/mol.
- Figure 5D/S13: UMA IRC paths for TSA–TSD and the activation-strain inputs below. [normalized_IRC_DFT_MLIP.csv](normalized_IRC_DFT_MLIP.csv) contains the tabulated DFT/UMA path comparison.
- Figure 5E/S20: [ts_bond_lengths.csv](ts_bond_lengths.csv) contains all-atom RMSDs and two forming-bond distances for TSA–TSD. The DFT TSB `r2` value is 2.060679 Å, rounding to 2.06 Å.

CSV labels `alkene`, `diazomethane`, `ts-a`–`ts-d`, and `prod-a`–`prod-d` correspond to S19's cyclopentenone, TMS-diazomethane, TSA–TSD, and Prod-A–D, respectively. The SI RMSDs for this study are all-atom values, distinct from the core heavy-atom benchmark.

## Starting geometries

The [starting geometries](starting_geometries/) use the same labels in their filenames and XYZ comments as the stationary-point CSVs. Their comments retain the recorded conformer energies, quasi-harmonic thermochemistry at 195.15 K and selection criterion. Transition-structure and product comments also retain their original conformer source labels. The 40-atom reactant is the alkene, C13H22O4Si.

| Starting XYZ | SI label |
| --- | --- |
| [alkene.xyz](starting_geometries/alkene.xyz) | Cyclopentenone / starting material 4 |
| [ts-a.xyz](starting_geometries/ts-a.xyz) | TSA |
| [prod-a.xyz](starting_geometries/prod-a.xyz) | Prod-A / product 5 |
| [ts-b.xyz](starting_geometries/ts-b.xyz) | TSB |
| [prod-b.xyz](starting_geometries/prod-b.xyz) | Prod-B / product 5b |
| [ts-c.xyz](starting_geometries/ts-c.xyz) | TSC |
| [prod-c.xyz](starting_geometries/prod-c.xyz) | Prod-C / product 5c |
| [ts-d.xyz](starting_geometries/ts-d.xyz) | TSD |
| [prod-d.xyz](starting_geometries/prod-d.xyz) | Prod-D / product 5d |

The four starting TS geometries provide the DFT comparison coordinates for Table S20 and Figure 5E. Their atom order matches the corresponding UMA frames in `cycloaddition_structures.xyz`.

## Activation-strain methods

For the TSA and TSB pathways, fragments were extracted from IRC geometries and their electronic energies evaluated at the same level as the ORCA UMA-S-1.1 jobs. The associated IRC tables record the model/task as `UMA-s-1p1/omol`. The decomposition uses electronic energies; Gibbs energies and thermal corrections are not used in the activation-strain terms.

The complex and both fragments are neutral singlets. [fragment_definition.csv](fragment_definition.csv) records the fragment membership using one-based indices in the full complex geometry. Fragment 1 has formula C13H22O4Si and fragment 2 has formula C4H10N2Si; their memberships together cover all 57 atoms of the complex without overlap.

The reaction coordinates are the two forming-bond distances, verified against the deposited IRC coordinates for both pathways:

| CSV coordinate | Atom pair, one-based | Units |
| --- | --- | --- |
| `r1_A` | C13–N14 | Å |
| `r2_A` | C2–C16 | Å |

At a sampled geometry `s`, the deposited electronic-energy decomposition is:

```text
strain_fragment1(s) = E_fragment1(s) - E_fragment1_reference
strain_fragment2(s) = E_fragment2(s) - E_fragment2_reference
strain_total(s)     = strain_fragment1(s) + strain_fragment2(s)
interaction(s)      = E_complex(s) - E_fragment1(s) - E_fragment2(s)
E_rel_separated(s)  = strain_total(s) + interaction(s)
```

The reference energies are those of the separated `alkene` and `diazomethane` entries in [gibbs_energies.csv](gibbs_energies.csv), specifically the electronic-energy column `E_hartree`:

| Fragment | Reference label | Electronic energy / hartree |
| --- | --- | --- |
| 1 | `alkene` | −1099.076497685 |
| 2 | `diazomethane` | −557.3972815208 |

Absolute electronic energies are stored in hartree. Strain, interaction, and relative energies are stored in kcal/mol; the table values reproduce the equations using 627.5094740631 kcal/mol per hartree. `closure_error_kcal` records the residual between the relative complex energy and the sum of the strain and interaction terms.

## Deposited path coverage

The activation-strain table contains 72 rows:

- **TSA:** 23 points, numbered 0–22. Their forming-bond distances match frames 23–1 of [TSA.xyz](irc/TSA.xyz), respectively, corresponding to IRC steps −1 through −23 in [TSA.csv](irc/TSA.csv). This covers the negative branch and excludes the step-0 transition state and positive branch.
- **TSB:** 49 points, numbered 0–48. Their forming-bond distances match frames 1–49 of [TSB.xyz](irc/TSB.xyz), respectively, corresponding to IRC steps −25 through +23 in [TSB.csv](irc/TSB.csv).

These mappings are established from the deposited bond distances. The CSV tables contain the complex and fragment electronic energies, thermochemical quantities, and activation-strain terms; the XYZ files provide the stationary-point, scan, and IRC geometries. Fragment membership is specified in `fragment_definition.csv` using one-based indices in the full complex geometry.
