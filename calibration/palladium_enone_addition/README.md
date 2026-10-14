# Pd conjugate-addition calibration — Figure 8B

This folder maps to manuscript Figure 8B, SI Tables S23/S24 and Figure S16. The 15 named entries match Table S3 (S-21–S-22) of the Q2MM study SI, DOI `10.1021/acs.joc.1c00136`. The reaction study is cited as DOI `10.1021/ja401713g`.

## Files and units

| File | Contents |
| --- | --- |
| [pd_series.csv](pd_series.csv) | 15 entries with structure labels, signed ee, experimental/raw predicted ΔΔG and experimental metadata fields |
| [pd_experimental_ee.csv](pd_experimental_ee.csv) | Signed experimental ee for P01–P15 in entry order |
| Thirty XYZ files | Competing `ts-r`/`ts-s` structures for P01–P15; coordinates in Å |
| [calibrate.py](calibrate.py) | Repeated-split linear fit using NumPy |
| [pd_train_test_flags.csv](pd_train_test_flags.csv) | 1,500 membership rows: 100 partitions ×15 entries |
| [pd_held_out_predictions.csv](pd_held_out_predictions.csv) | 700 held-out predictions: 100 partitions ×7 test entries |
| [pd_split_metrics.csv](pd_split_metrics.csv) | 100 rows of memberships, coefficients and raw/calibrated test RMSE/MAE; supplies the paired MAE arrays for S16 |
| [calibration_summary.csv](calibration_summary.csv) | Coefficient/error means and sample SDs, pooled held-out errors, unique-entry raw errors and full-data fit |

ΔΔG, residuals, intercepts, RMSEs and MAEs are in kcal/mol; slope is dimensionless. Temperature fields use kelvin. Across-partition uncertainties are sample SDs (`ddof=1`). Split-average, pooled held-out and unique-entry statistics have distinct definitions.

## Running the calibration

Keep `calibrate.py` beside `pd_series.csv` and run:

```sh
python -m pip install numpy
python calibrate.py
```

The script and saved analyses use **100 unique partitions**, eight training/seven test entries and seed `20260810`. Training combinations are sampled without replacement. Row order determines the seeded selection, so retain the supplied input order. Inputs and outputs are resolved beside the script, including when launched from another working directory. A rerun overwrites the derived output CSVs and writes `runtime.json` with the NumPy version and split settings.

The fit is `ddG_exp = slope × ddG_uma_raw + intercept`, using training entries only. Raw and calibrated test errors use the held-out entries. There is no single train/test assignment for a plot aggregated over repeated partitions.

The saved split-average calibrated RMSE is **0.670552 ±0.132768 kcal/mol** and MAE **0.566282 ±0.139118 kcal/mol**. Calibrated MAE is lower than raw MAE in **86/100 partitions**; the reduction in mean held-out MAE is **27.0%**.

P03 (`pyox1-11-12`) corresponds to `enone031`, −87% ee; P11 (`pyox1-12-1`) corresponds to `enone047`, +87% ee. These values are recorded consistently in the series and experimental-ee tables.
