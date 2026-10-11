# Cycloaddition reaction landscape

This directory contains the [ten stationary-point geometries](cycloaddition_structures.xyz), [electronic and Gibbs energies](gibbs_energies.csv), [two-dimensional scan](2d-scan/scan2d.csv), [IRC trajectories and energy tables](irc/), and [activation-strain decomposition](activation_strain.csv).

## Activation-strain methods

For the TSA and TSB pathways, fragments were extracted from IRC geometries and their electronic energies evaluated at the same level as the ORCA UMA-S-1.1 jobs, as confirmed by the repository author. The associated IRC tables record the model/task as `UMA-s-1p1/omol`. The decomposition uses electronic energies; Gibbs energies and thermal corrections are not used in the activation-strain terms.

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

These mappings are established from the deposited bond distances. The reason for the different path coverage is not yet documented.

## Remaining provenance

The table includes source output filenames and line numbers, but the named source outputs and fragment manifests are not deposited. The precise ORCA version, external-model setup, calculation inputs, fragment relaxation/constraint treatment, and reference-fragment optimization protocol remain to be recorded. No unreported settings or path-selection rationale have been inferred.
