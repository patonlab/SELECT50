Standalone NI calibration

Move this directory anywhere; keep calibrate.py and ni_series.csv together.
Install: python -m pip install numpy
Run: python calibrate.py
Only NumPy is required. Inputs and outputs are resolved beside the script, so it
also works when launched from another working directory. Existing derived outputs
are overwritten on rerun; the input series CSV is not changed.

Default seed: 20260810. 100 unique splits, 7 training and 6 test
entries. Row order determines seed reproducibility; retain the supplied order.
The constants at the top of calibrate.py control these settings. Both deposited
Ni and Pd output tables contain 100 splits. The Pd script currently requests
1000, while the manuscript/SI say 100; the original figure protocol is unresolved.
Changing sample size changes the selected partitions. See the Pd README.md.

Fit: ddG_exp = slope * ddG_uma_raw + intercept, using training entries only.
Held-out errors use test entries only. Raw errors compare ddG_uma_raw directly
with ddG_exp. All energies and errors are kcal/mol; slope is dimensionless.
Slope/intercept and splitwise error SD use sample SD across splits (ddof=1).
There is no single train/test flag for an aggregate repeated-split figure.

Outputs:
- ni_train_test_flags.csv: train/test assignment for every entry in every split.
- ni_split_metrics.csv: memberships, fitted coefficients, raw/calibrated test errors.
- ni_held_out_predictions.csv: raw/calibrated held-out predictions and residuals.
- calibration_summary.csv: means, SDs, pooled errors, unique-entry raw errors, full-data fit.
- runtime.json: NumPy version, seed, training size and number of splits.

Input columns: entry_id, ligand_or_substrate, ee_percent, temperature_K,
ddG_exp, ddG_uma_raw, reference. Ligand identifiers were taken from the deposited
geometry filenames. Original entry order and energy values are preserved.

Experimental metadata gaps in the input remain blank. The separately supplied
working SI Table S24 contains all 13 Ni ee values in entry order:
59, 90, 85, 65, 66, 89, 87, 63, 64, 85, 84, 75, 86 percent.
The corresponding geometry ligand labels are l1, l4, l7, l9, l11, l12, l14, l20,
l23, l26, l27, l28, l29. These values are available for the later experimental
metadata task; they have not been inserted into ni_series.csv here. The S24
row-to-ligand mapping still needs source confirmation.

The drafts cite main reference 90/SI reference 87, DOI 10.1021/jacs.1c08105.
Per-entry literature attribution and actual experimental temperatures remain
unverified. Agreement of ee-derived dG at 298.15 K with rounded CSV values does
not establish the experimental temperature.

Ni P12 retains its original raw calibration dG (3.727905467 kcal/mol). Its XYZ
comments imply an electronic-energy separation of about 1.768 kcal/mol; these
are different energy types, so confirm the Gibbs correction/source before
treating the difference as an error. Original thermochemistry records are absent.

The saved arithmetic has been checked using its memberships; the script does
not generate chemistry calculations, plots, or all printed R2/Spearman metrics.
runtime.json is written on rerun but is not currently deposited.
