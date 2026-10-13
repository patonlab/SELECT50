#!/usr/bin/env python3
"""Standalone PD repeated-split OLS calibration and raw baselines.
Run: python calibrate.py
Only NumPy is required. Input energy differences and output errors are kcal/mol.
"""
from pathlib import Path
import csv
import itertools
import json
import numpy as np

ROOT = Path(__file__).absolute().parent
SERIES = 'pd'
RANDOM_SEED = 20260810
EXPECTED_ROWS = 15
TRAIN_SIZE = 8
N_SPLITS = 100


def read_csv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def errors(observed, predicted):
    residual = observed - predicted
    return float(np.sqrt(np.mean(residual ** 2))), float(np.mean(np.abs(residual)))


def run_series(series, config, seed):
    rows = read_csv(ROOT / f'{series}_series.csv')
    ids = [row['entry_id'] for row in rows]
    x = np.array([float(row['ddG_uma_raw']) for row in rows])
    y = np.array([float(row['ddG_exp']) for row in rows])
    if len(rows) != config['n_rows'] or len(set(ids)) != len(ids):
        raise ValueError(f'{series}: unexpected row count or duplicate IDs')
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError(f'{series}: energy differences must be finite')
    if not (2 <= config['train_size'] <= len(rows) - 2):
        raise ValueError('Leave at least two entries in each training and test set')
    combinations = list(itertools.combinations(range(len(rows)), config['train_size']))
    if not (2 <= config['n_splits'] <= len(combinations)):
        raise ValueError('Split count must be between 2 and the number of unique partitions')
    # Preserve the original choice without replacement over lexicographic combinations.
    choices = np.random.default_rng(seed).choice(
        len(combinations), size=config['n_splits'], replace=False)
    split_rows, flags, held_out = [], [], []
    for number, choice in enumerate(choices):
        train = np.array(combinations[int(choice)])
        mask = np.zeros(len(rows), dtype=bool)
        mask[train] = True
        test = np.flatnonzero(~mask)
        if np.ptp(x[train]) == 0:
            raise ValueError(f'{series}: training raw energies have no variation')
        slope, intercept = np.polyfit(x[train], y[train], 1)
        prediction = slope * x[test] + intercept
        raw_rmse, raw_mae = errors(y[test], x[test])
        calibrated_rmse, calibrated_mae = errors(y[test], prediction)
        split = f'split_{number + 1:04d}'
        train_ids = ', '.join(ids[i] for i in train)
        test_ids = ', '.join(ids[i] for i in test)
        split_rows.append(dict(
            split=split, train_IDs=train_ids, test_IDs=test_ids,
            slope=float(slope), intercept=float(intercept),
            raw_test_RMSE=raw_rmse, raw_test_MAE=raw_mae,
            calibrated_test_RMSE=calibrated_rmse, calibrated_test_MAE=calibrated_mae))
        for i, entry_id in enumerate(ids):
            flags.append(dict(split=split, entry_id=entry_id,
                              train_test='train' if mask[i] else 'test'))
        for i, predicted in zip(test, prediction):
            held_out.append(dict(split=split, entry_id=ids[i], ddG_exp=float(y[i]),
                                 ddG_uma_raw=float(x[i]),
                                 ddG_calibrated=float(predicted),
                                 raw_residual_exp_minus_prediction=float(y[i]-x[i]),
                                 calibrated_residual_exp_minus_prediction=float(y[i]-predicted)))
    write_csv(ROOT / f'{series}_split_metrics.csv', split_rows)
    write_csv(ROOT / f'{series}_train_test_flags.csv', flags)
    write_csv(ROOT / f'{series}_held_out_predictions.csv', held_out)
    summary = dict(series=series, seed=seed, n_entries=len(rows),
                   n_splits=len(split_rows), train_size=config['train_size'],
                   test_size=len(rows)-config['train_size'])
    # SD is sample SD across training-set fits/split metrics, not a standard error.
    for key in ['slope', 'intercept', 'raw_test_RMSE', 'raw_test_MAE',
                'calibrated_test_RMSE', 'calibrated_test_MAE']:
        values = np.array([row[key] for row in split_rows])
        summary[key+'_mean'] = float(np.mean(values))
        summary[key+'_sd'] = float(np.std(values, ddof=1))
    observed = np.array([row['ddG_exp'] for row in held_out])
    raw = np.array([row['ddG_uma_raw'] for row in held_out])
    calibrated = np.array([row['ddG_calibrated'] for row in held_out])
    summary['raw_pooled_held_out_RMSE'], summary['raw_pooled_held_out_MAE'] = errors(observed, raw)
    summary['calibrated_pooled_held_out_RMSE'], summary['calibrated_pooled_held_out_MAE'] = errors(observed, calibrated)
    # One raw prediction per entry: removes unequal repeated-test exposure weighting.
    summary['raw_unique_entry_RMSE'], summary['raw_unique_entry_MAE'] = errors(y, x)
    slope, intercept = np.polyfit(x, y, 1)
    summary['full_data_slope'], summary['full_data_intercept'] = float(slope), float(intercept)
    return summary


def main():
    config = dict(n_rows=EXPECTED_ROWS, train_size=TRAIN_SIZE, n_splits=N_SPLITS)
    summaries = [run_series(SERIES, config, RANDOM_SEED)]
    write_csv(ROOT / 'calibration_summary.csv', summaries)
    (ROOT / 'runtime.json').write_text(json.dumps({'numpy_version': np.__version__, 'seed': RANDOM_SEED,
                                                 'train_size': TRAIN_SIZE, 'n_splits': N_SPLITS}, indent=2)+'\n')
    for row in summaries:
        print(f"{row['series'].upper()}: {row['n_splits']} splits, "
              f"slope {row['slope_mean']:.6f} +/- {row['slope_sd']:.6f}, "
              f"intercept {row['intercept_mean']:.6f} +/- {row['intercept_sd']:.6f}; "
              f"raw test RMSE {row['raw_test_RMSE_mean']:.6f} +/- {row['raw_test_RMSE_sd']:.6f}, "
              f"MAE {row['raw_test_MAE_mean']:.6f} +/- {row['raw_test_MAE_sd']:.6f}")
    print(f'Wrote split memberships, predictions and calibration summary to {ROOT}')


if __name__ == '__main__':
    main()
