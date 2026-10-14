Ni epoxide-arylation calibration — Figure 8A

This dataset maps to manuscript Figure 8A, SI Tables S23/S24 and Figure S15.
The reaction-study DOI is 10.1021/jacs.1c08105.

The author-confirmed experimental temperature is 34 °C (307.15 K).
Experimental gaps are converted from signed ee using
ddG_exp = R*T*ln((1+ee/100)/(1-ee/100)), with
R = 0.00198720425864083 kcal mol^-1 K^-1. Raw UMA gaps are retained.

Move this directory anywhere; keep calibrate.py and ni_series.csv together.
Install: python -m pip install numpy
Run: python calibrate.py
Only NumPy is required. Inputs and outputs are resolved beside the script,
including when launched from another working directory. Existing derived outputs
are overwritten on rerun; the input series CSV is retained.

Default seed: 20260810. The script and saved analysis use 100 unique splits,
with 7 training and 6 test entries from a 13-entry series. Training combinations
are sampled without replacement. Row order determines seed reproducibility;
retain the supplied order. Both deposited Ni and Pd analyses use 100 partitions.
The constants at the top of calibrate.py control these settings.

Fit: ddG_exp = slope * ddG_uma_raw + intercept, using training entries only.
Held-out errors use test entries only. Raw errors compare ddG_uma_raw directly
with ddG_exp. All energies and errors are kcal/mol; slope is dimensionless.
Slope/intercept and splitwise error SD use sample SD across splits (ddof=1).
There is no single train/test flag for an aggregate repeated-split figure.

Files:
- ni_series.csv: entry labels, experimental/raw UMA gaps and metadata fields.
- ni_experimental_ee.csv: signed experimental ee for the 13 entries in order.
- ni_train_test_flags.csv: train/test assignments for every entry in every split.
- ni_split_metrics.csv: memberships, coefficients and raw/calibrated test errors.
- ni_held_out_predictions.csv: raw/calibrated held-out predictions and residuals.
- calibration_summary.csv: means, SDs, pooled errors, unique-entry raw errors
  and full-data fit.
- ni_calibration_plot_data.csv: 100 paired raw/calibrated held-out MAEs for
  Figure S15, test R2/Spearman values, and mean-prediction null errors.
- ni_calibration_plot_points.csv: 13 experimental gaps, per-entry held-out
  prediction means/sample SDs and test counts for Figure 8A.
- calibration_refit.json: conversion temperature/constant, input and membership
  hashes, retained seed and refit implementation.
- runtime.json: written on rerun; records NumPy version, seed, training size
  and number of splits.

Input columns: entry_id, ligand_or_substrate, ee_percent, temperature_K,
ddG_exp, ddG_uma_raw, reference. Ligand identifiers in the series correspond
to the deposited geometry filenames.

The experimental-ee table contains, in entry order:
59, 90, 85, 65, 66, 89, 87, 63, 64, 85, 84, 75, 86 percent.
The series geometry labels are l1, l4, l7, l9, l11, l12, l14, l20,
l23, l26, l27, l28 and l29.

The scripts fit the supplied Gibbs-energy differences; XYZ-comment electronic
energies are a separate quantity. The script generates calibration tables
and a runtime record.

To refresh the experimental conversion and all plot/statistic tables while
retaining the saved 100 memberships, run from the repository root:
python calibration/refresh_saved_calibration.py
This shared driver uses only the Python standard library and performs no new
random sampling. It also updates repeated_random_split_metrics.csv; that table's
test error bands are ±0.10 and ±0.25 kcal/mol, as displayed in Figure 8A.

At 307.15 K, calibrated test RMSE is 0.148102 ±0.028717 kcal/mol and MAE is
0.123077 ±0.030579 kcal/mol. MAE improves in all 100 partitions; the reduction
in mean held-out MAE is 93.6%.

For Figure S15, plot raw_test_MAE on x and calibrated_test_MAE on y. The
training_mean_null_test_* columns predict the mean of the training experimental
gaps. The test_mean_null_* columns use the held-out mean and describe the
benchmark underlying R2 = 1 - SSE/SST. Raw UMA is a separate baseline.
Against the training-mean null, calibrated MAE improves in 100/100 partitions;
the reduction in mean held-out MAE is 64.9%.
