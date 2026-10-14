# Palladium conjugate-addition calibration

This dataset contains 15 entries and corresponds to Figure 8B, SI Tables S23–S24 and Figure S16. The 15 entries correspond to Table S3 (S-21–S-22) of the Q2MM study SI, DOI `10.1021/acs.joc.1c00136`. The reaction-study DOI is `10.1021/ja401713g`.

## Files

| File | Contents |
| --- | --- |
| [pd_series.csv](pd_series.csv) | Entry IDs, structure labels, signed ee, experimental and UMA-S-1.2 Gibbs gaps, temperature and references |
| [pd_experimental_ee.csv](pd_experimental_ee.csv) | Signed experimental ee values in entry order |
| [pd_train_test_flags.csv](pd_train_test_flags.csv) | Training/test assignments for all 15 entries in each of 100 partitions |
| [pd_split_metrics.csv](pd_split_metrics.csv) | Memberships, fitted coefficients and raw/calibrated held-out errors |
| [pd_held_out_predictions.csv](pd_held_out_predictions.csv) | 700 held-out predictions and residuals |
| [calibration_summary.csv](calibration_summary.csv) | Mean and sample SD of partition metrics, pooled errors, null comparisons and full-data fit |
| [pd_calibration_plot_data.csv](pd_calibration_plot_data.csv) | 100 paired MAEs, test R², Spearman correlation and mean-prediction null errors for Figure S16 |
| [pd_calibration_plot_points.csv](pd_calibration_plot_points.csv) | Per-entry held-out prediction means, sample SDs and test counts for Figure 8B |
| [calibration_refit.json](calibration_refit.json) | Conversion conditions, partition seed, fitting method and input hashes |

The XYZ files provide the competing transition structures identified in the series table. Coordinates use ångströms. Gibbs gaps, residuals, intercepts, RMSE and MAE use kcal/mol; slopes are dimensionless.

## Calibration

Experimental gaps are converted from signed ee at **333 K**:

```text
ddG_exp = R*T*ln((1+ee/100)/(1-ee/100))
R = 0.00198720425864083 kcal mol^-1 K^-1
```

Each of the 100 recorded partitions has 8 training and 7 test entries. The partition-selection seed is `20260810`. The model `ddG_exp = slope*ddG_uma_raw + intercept` is fitted on training entries and evaluated on held-out entries. Across-partition and per-entry uncertainties are sample standard deviations.

The predictive null uses the mean experimental gap of the training set. The separate test-mean baseline defines `R² = 1 - SSE/SST`. Raw UMA predictions provide a third comparison. For Figure S16, plot `raw_test_MAE` against `calibrated_test_MAE`; Figure 8B uses the per-entry plotting table.

## Run

From the repository root:

```sh
python3 calibration/palladium_enone_addition/calibrate.py
```

The command uses only the Python standard library and regenerates all derived CSVs and calibration metadata for this series using its recorded partitions. Keep the repository's calibration directory structure intact. To analyze both series together, run:

```sh
python3 calibration/refresh_saved_calibration.py
```
