# Pd conjugate-addition calibration — Figure 9B

This folder maps to manuscript Figure 9B and SI Tables S23/S24 and Figure S16 in the separately supplied working drafts. Literature attribution in those drafts is main reference 91/SI reference 88, DOI `10.1021/ja401713g`. Per-entry experimental attribution has not been verified against the publication.

## Deposited files and units

| File | Contents |
| --- | --- |
| [pd_series.csv](pd_series.csv) | 15 input entries: ligand/substrate filename label, signed ee percent, experimental/raw UMA ΔΔG, incomplete temperature/reference fields |
| Thirty XYZ files | Competing `ts-r`/`ts-s` structures for P01–P15; coordinates in Å |
| [calibrate.py](calibrate.py) | NumPy repeated-split linear fit, experimental ΔΔG = slope × raw predicted ΔΔG + intercept, fitted on training entries only |
| [pd_train_test_flags.csv](pd_train_test_flags.csv) | 1,500 membership rows: 100 partitions ×15 entries |
| [pd_held_out_predictions.csv](pd_held_out_predictions.csv) | 700 held-out predictions: 100 partitions ×7 test entries |
| [pd_split_metrics.csv](pd_split_metrics.csv) | 100 rows of memberships, coefficients, and raw/calibrated test RMSE/MAE; supplies paired MAE arrays for S16 |
| [calibration_summary.csv](calibration_summary.csv) | Saved coefficient/error means and sample SDs, pooled held-out errors, unique-entry raw errors, and full-data fit |

ΔΔG, residuals, intercepts, RMSEs, and MAEs are in kcal/mol; slope is dimensionless. Temperature, when present, is in kelvin. A sample SD across partitions is distinct from a standard error. Split-average, pooled held-out, and unique-entry statistics are distinct summary definitions.

## Reproduction discrepancy

The saved analysis has **100 unique partitions**, eight training/seven test entries, and seed `20260810`. The current script instead sets **`N_SPLITS = 1000`**. It chooses unique index combinations without replacement; requested sample size affects selected memberships. Running it now overwrites the deposited output files and produces a different partition sample. It also writes `runtime.json`, which is not currently deposited.

The script requires NumPy and resolves files beside itself; retain it with `pd_series.csv`. Resolve the final split count and experimental records before regenerating saved results. The saved 100-split summaries give calibrated RMSE **0.670552 ±0.132768** and MAE **0.566282 ±0.139118 kcal/mol**. The working manuscript/Figure 9/S23 instead report **0.64 ±0.12** and **0.53 ±0.12**. The saved arithmetic is internally consistent; its agreement with the publication's selected input/protocol remains unresolved.

## Experimental and figure gaps

- All experimental temperatures and references are blank in the input CSV. Its ee/ΔΔG pairs algebraically match conversion at 298.15 K; that is not evidence of the original experimental temperature.
- P03 (`pyox1-11-12`) is −87% ee in the CSV versus −92% in S24. P11 (`pyox1-12-1`) is +87% versus +81%. Keep these conflicts visible until the original experimental entries are identified; corrections must propagate to ΔΔG and derived statistics.
- Figure 9's caption describes blue training circles/red held-out squares, while its current image shows red circular predictions with error bars. Define whether points aggregate repeated held-out predictions or represent one split. An aggregate plot has no single train/test flag.
- Current scripts do not generate the plots or every R²/Spearman statistic shown in Figure 9/S23. Original plotting/metric provenance is still required.
- S24 contains the chemical structures as ChemDraw OLE objects, with incomplete identity fields. Current S13/S14 are IRC figures; S15/S16 are the Ni/Pd improvement plots. Main-text references to S17/S18 are stale.
- The actual UMA checkpoint and raw predicted Gibbs-energy/thermochemistry provenance remain to be documented. No input energy or experimental value was changed during the documentation update.
