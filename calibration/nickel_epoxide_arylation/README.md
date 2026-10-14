# Nickel epoxide-arylation calibration

This dataset contains 13 entries and corresponds to Figure 8A, SI Tables S23–S24 and Figure S15. The reaction-study DOI is `10.1021/jacs.1c08105`.

## Files

| File | Contents |
| --- | --- |
| [ni_series.csv](ni_series.csv) | Entry IDs, structure labels, signed ee, experimental and UMA-S-1.2 Gibbs gaps, temperature and references |
| [ni_experimental_ee.csv](ni_experimental_ee.csv) | Signed experimental ee values in entry order |
| [ni_train_test_flags.csv](ni_train_test_flags.csv) | Training/test assignments for all 13 entries in each of 100 partitions |
| [ni_split_metrics.csv](ni_split_metrics.csv) | Memberships, fitted coefficients and raw/calibrated held-out errors |
| [ni_held_out_predictions.csv](ni_held_out_predictions.csv) | 600 held-out predictions and residuals |
| [calibration_summary.csv](calibration_summary.csv) | Mean and sample SD of partition metrics, pooled errors, null comparisons and full-data fit |
| [ni_calibration_plot_data.csv](ni_calibration_plot_data.csv) | 100 paired MAEs, test R², Spearman correlation and mean-prediction null errors for Figure S15 |
| [ni_calibration_plot_points.csv](ni_calibration_plot_points.csv) | Per-entry held-out prediction means, sample SDs and test counts for Figure 8A |
| [repeated_random_split_metrics.csv](repeated_random_split_metrics.csv) | Detailed training/test metrics, pairwise ordering and ±0.10/±0.25 kcal/mol error bands |
| [calibration_refit.json](calibration_refit.json) | Conversion conditions, partition seed, fitting method and input hashes |

The XYZ files provide the competing transition structures identified in the series table. Coordinates use ångströms. Gibbs gaps, residuals, intercepts, RMSE and MAE use kcal/mol; slopes are dimensionless.

## Calibration

Experimental gaps are converted from signed ee at **307.15 K (34 °C)**:

```text
ddG_exp = R*T*ln((1+ee/100)/(1-ee/100))
R = 0.00198720425864083 kcal mol^-1 K^-1
```

Each of the 100 recorded partitions has 7 training and 6 test entries. The partition-selection seed is `20260810`. The model `ddG_exp = slope*ddG_uma_raw + intercept` is fitted on training entries and evaluated on held-out entries. Across-partition and per-entry uncertainties are sample standard deviations.

The predictive null uses the mean experimental gap of the training set. The separate test-mean baseline defines `R² = 1 - SSE/SST`. Raw UMA predictions provide a third comparison. For Figure S15, plot `raw_test_MAE` against `calibrated_test_MAE`; Figure 8A uses the per-entry plotting table.

## Run

From the repository root:

```sh
python3 calibration/nickel_epoxide_arylation/calibrate.py
```

The command uses only the Python standard library and regenerates all derived CSVs and calibration metadata for this series using its recorded partitions. Keep the repository's calibration directory structure intact. To analyze both series together, run:

```sh
python3 calibration/refresh_saved_calibration.py
```
