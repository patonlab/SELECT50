"""Analyze Ni/Pd calibration using the recorded training/test partitions.

Run from the repository root: python calibration/refresh_saved_calibration.py
Only the Python standard library is required. No random sampling is performed.
"""
from pathlib import Path
import collections
import csv
import hashlib
import io
import itertools
import json
import math
import platform
import statistics as st

ROOT = Path(__file__).resolve().parent
R = 0.00198720425864083  # kcal mol^-1 K^-1
FAMILIES = {'ni': ('nickel_epoxide_arylation', 13, 7),
            'pd': ('palladium_enone_addition', 15, 8)}


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def ee_to_ddg(ee_percent, temperature):
    if not math.isfinite(ee_percent) or not -100 < ee_percent < 100:
        raise ValueError('Signed ee must be finite and strictly between -100 and 100')
    if not math.isfinite(temperature) or temperature <= 0:
        raise ValueError('Temperature must be finite and positive in kelvin')
    ee = ee_percent / 100
    return R * temperature * math.log((1 + ee) / (1 - ee))


def fit(x, y):
    xm, ym = st.mean(x), st.mean(y)
    denominator = sum((v - xm) ** 2 for v in x)
    if denominator == 0:
        raise ValueError('Training predictions have no variation')
    slope = sum((a - xm) * (b - ym) for a, b in zip(x, y)) / denominator
    return slope, ym - slope * xm


def ranks(values):
    groups = collections.defaultdict(list)
    for rank, index in enumerate(sorted(range(len(values)), key=values.__getitem__), 1):
        groups[values[index]].append((rank, index))
    result = [0.] * len(values)
    for group in groups.values():
        average = st.mean(rank for rank, _ in group)
        for _, index in group:
            result[index] = average
    return result


def scores(observed, predicted):
    errors = [a - b for a, b in zip(observed, predicted)]
    sst = sum((v - st.mean(observed)) ** 2 for v in observed)
    if not sst:
        raise ValueError('Observed test values have no variation')
    x, y = ranks(observed), ranks(predicted)
    xm, ym = st.mean(x), st.mean(y)
    denominator = math.sqrt(sum((v - xm) ** 2 for v in x) * sum((v - ym) ** 2 for v in y))
    rho = sum((a - xm) * (b - ym) for a, b in zip(x, y)) / denominator if denominator else None
    return {'RMSE': math.sqrt(st.mean(e * e for e in errors)),
            'MAE': st.mean(abs(e) for e in errors),
            'R2': 1 - sum(e * e for e in errors) / sst, 'rho': rho}


def table_bytes(path, rows):
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    bom = b'\xef\xbb\xbf' if path.exists() and path.read_bytes().startswith(b'\xef\xbb\xbf') else b''
    return bom + stream.getvalue().encode('utf-8')


def comparable_pairs(observed, predicted):
    pairs = [(observed[i] - observed[j], predicted[i] - predicted[j])
             for i, j in itertools.combinations(range(len(observed)), 2)
             if observed[i] != observed[j]]
    correct = sum(a * b > 0 for a, b in pairs)
    return correct, len(pairs), correct / len(pairs)


def prepare(family, folder, n, ntrain):
    base = ROOT / folder
    series_path = base / f'{family}_series.csv'
    series = read(series_path)
    ids = [r['entry_id'] for r in series]
    if len(series) != n or len(set(ids)) != n:
        raise ValueError(f'{family}: invalid entry coverage')
    ee = read(base / f'{family}_experimental_ee.csv')
    if [int(r['ligand_id']) for r in ee] != list(range(1, n + 1)):
        raise ValueError(f'{family}: unexpected experimental-ee ordering')
    for row, source in zip(series, ee):
        if row['entry_id'] != f"P{int(source['ligand_id']):02d}" or float(row['ee_percent']) != float(source['experimental_ee_percent_signed']):
            raise ValueError(f'{family}: ee/entry mapping differs')
        row['ddG_exp'] = str(ee_to_ddg(float(row['ee_percent']), float(row['temperature_K'])))
    data = {r['entry_id']: r for r in series}
    if any(not math.isfinite(float(r['ddG_uma_raw'])) for r in series):
        raise ValueError(f'{family}: non-finite raw prediction')
    old_splits = read(base / f'{family}_split_metrics.csv')
    if len(old_splits) != 100 or len({r['split'] for r in old_splits}) != 100:
        raise ValueError(f'{family}: expected 100 saved partitions')
    split_rows, held, plot, train_scores, test_details, membership = [], [], [], {}, {}, set()
    for row in old_splits:
        train = [v.strip() for v in row['train_IDs'].split(',')]
        test = [v.strip() for v in row['test_IDs'].split(',')]
        if len(train) != ntrain or len(test) != n - ntrain or set(train) & set(test) or set(train + test) != set(ids) or tuple(train) in membership:
            raise ValueError(f"{family}/{row['split']}: invalid membership")
        membership.add(tuple(train))
        xtrain = [float(data[i]['ddG_uma_raw']) for i in train]
        ytrain = [float(data[i]['ddG_exp']) for i in train]
        slope, intercept = fit(xtrain, ytrain)
        x = [float(data[i]['ddG_uma_raw']) for i in test]
        y = [float(data[i]['ddG_exp']) for i in test]
        prediction = [slope * v + intercept for v in x]
        raw, cal = scores(y, x), scores(y, prediction)
        training_null = scores(y, [st.mean(ytrain)] * len(test))
        test_null = scores(y, [st.mean(y)] * len(test))
        split = row['split']
        split_rows.append(dict(split=split, train_IDs=row['train_IDs'], test_IDs=row['test_IDs'],
                               slope=slope, intercept=intercept,
                               raw_test_RMSE=raw['RMSE'], raw_test_MAE=raw['MAE'],
                               calibrated_test_RMSE=cal['RMSE'], calibrated_test_MAE=cal['MAE']))
        plot.append(dict(split=split, raw_test_MAE=raw['MAE'], calibrated_test_MAE=cal['MAE'],
                         improved_vs_raw=int(cal['MAE'] < raw['MAE']),
                         training_mean_null_test_MAE=training_null['MAE'],
                         training_mean_null_test_RMSE=training_null['RMSE'],
                         test_mean_null_MAE=test_null['MAE'], test_mean_null_RMSE=test_null['RMSE'],
                         test_R2=cal['R2'], test_spearman_rho=cal['rho']))
        train_scores[split] = (ytrain, [slope * v + intercept for v in xtrain])
        test_details[split] = (y, prediction)
        for entry, yy, xx, pp in zip(test, y, x, prediction):
            held.append(dict(split=split, entry_id=entry, ddG_exp=yy, ddG_uma_raw=xx, ddG_calibrated=pp,
                             raw_residual_exp_minus_prediction=yy - xx,
                             calibrated_residual_exp_minus_prediction=yy - pp))
    expected_flags = {(row['split'], i): ('train' if i in row['train_IDs'].split(', ') else 'test')
                      for row in old_splits for i in ids}
    flags = read(base / f'{family}_train_test_flags.csv')
    if len(flags) != n * 100 or {(r['split'], r['entry_id']): r['train_test'] for r in flags} != expected_flags:
        raise ValueError(f'{family}: saved membership flags differ')
    summary = dict(read(base / 'calibration_summary.csv')[0])
    for key in ('slope', 'intercept', 'raw_test_RMSE', 'raw_test_MAE', 'calibrated_test_RMSE', 'calibrated_test_MAE'):
        values = [r[key] for r in split_rows]
        summary[key + '_mean'], summary[key + '_sd'] = st.mean(values), st.stdev(values)
    for label, field in [('raw', 'ddG_uma_raw'), ('calibrated', 'ddG_calibrated')]:
        pooled = scores([r['ddG_exp'] for r in held], [r[field] for r in held])
        for metric in ('RMSE', 'MAE'):
            summary[f'{label}_pooled_held_out_{metric}'] = pooled[metric]
    x = [float(r['ddG_uma_raw']) for r in series]
    y = [float(r['ddG_exp']) for r in series]
    raw_unique = scores(y, x)
    summary['raw_unique_entry_RMSE'], summary['raw_unique_entry_MAE'] = raw_unique['RMSE'], raw_unique['MAE']
    summary['full_data_slope'], summary['full_data_intercept'] = fit(x, y)
    for key in ('test_R2', 'test_spearman_rho', 'training_mean_null_test_MAE', 'training_mean_null_test_RMSE'):
        values = [r[key] for r in plot]
        summary[key + '_mean'], summary[key + '_sd'] = st.mean(values), st.stdev(values)
    summary['n_improved_vs_raw'] = sum(r['improved_vs_raw'] for r in plot)
    summary['MAE_reduction_percent'] = 100 * (1 - summary['calibrated_test_MAE_mean'] / summary['raw_test_MAE_mean'])
    summary['n_improved_vs_training_mean_null'] = sum(
        r['calibrated_test_MAE'] < r['training_mean_null_test_MAE'] for r in plot)
    summary['MAE_reduction_vs_training_mean_null_percent'] = 100 * (
        1 - summary['calibrated_test_MAE_mean'] / summary['training_mean_null_test_MAE_mean'])
    summary['temperature_K'] = float(series[0]['temperature_K'])
    if len({r['temperature_K'] for r in series}) != 1:
        raise ValueError(f'{family}: expected a common experimental temperature')
    by_entry = collections.defaultdict(list)
    for row in held:
        by_entry[row['entry_id']].append(row['ddG_calibrated'])
    points = [dict(entry_id=r['entry_id'], ligand_or_substrate=r['ligand_or_substrate'],
                   ddG_exp=r['ddG_exp'], ddG_uma_raw=r['ddG_uma_raw'],
                   calibrated_prediction_mean=st.mean(by_entry[r['entry_id']]),
                   calibrated_prediction_sd=st.stdev(by_entry[r['entry_id']]),
                   n_held_out_partitions=len(by_entry[r['entry_id']])) for r in series]
    outputs = {series_path: series, base / f'{family}_split_metrics.csv': split_rows,
               base / f'{family}_held_out_predictions.csv': held,
               base / 'calibration_summary.csv': [summary],
               base / f'{family}_calibration_plot_data.csv': plot,
               base / f'{family}_calibration_plot_points.csv': points}
    extra = base / 'repeated_random_split_metrics.csv'
    if extra.exists():
        detailed = read(extra)
        splits = {r['split']: r for r in split_rows}
        if len(detailed) != 100 or {r['split'] for r in detailed} != set(splits):
            raise ValueError('Unexpected detailed Ni partition coverage')
        for row in detailed:
            split = row['split']
            if (row['train_IDs'], row['test_IDs']) != (splits[split]['train_IDs'], splits[split]['test_IDs']):
                raise ValueError('Detailed Ni memberships differ')
            row.update(slope=splits[split]['slope'], intercept=splits[split]['intercept'])
            for label, pair in [('train', train_scores[split]), ('test', test_details[split])]:
                observed, predicted = pair
                metric = scores(observed, predicted)
                row.update({label + '_R2': metric['R2'], label + '_RMSE_kcal_mol': metric['RMSE'],
                            label + '_MAE_kcal_mol': metric['MAE'], label + '_spearman_rho': metric['rho']})
                correct, comparable, accuracy = comparable_pairs(observed, predicted)
                row[label + '_pairwise_ordering_accuracy'] = accuracy
                if label == 'test':
                    errors = [abs(a - b) for a, b in zip(observed, predicted)]
                    row.update(test_max_abs_error_kcal_mol=max(errors),
                               test_n_correctly_ordered_pairs=correct, test_n_comparable_pairs=comparable,
                               test_n_within_inner=sum(e <= .10 for e in errors),
                               test_n_within_outer=sum(e <= .25 for e in errors))
                    # For six test entries, the usual two-sided t approximation has df=4.
                    rho = abs(metric['rho'])
                    row['test_spearman_p_value'] = (1 - rho) ** 2 * (1 + rho / 2)
        outputs[extra] = detailed
    provenance = dict(temperature_K=summary['temperature_K'],
                      temperature_source=('Reaction study, DOI 10.1021/jacs.1c08105' if family == 'ni'
                                          else 'Reaction study, DOI 10.1021/ja401713g'),
                      R_kcal_mol_K=R, formula='R*T*ln((1+ee/100)/(1-ee/100))',
                      seed=int(summary['seed']), n_splits=100,
                      partition_selection='100 recorded training/test partitions',
                      input_sha256=hashlib.sha256(table_bytes(series_path, series)).hexdigest(),
                      membership_flags_sha256=hashlib.sha256((base / f'{family}_train_test_flags.csv').read_bytes()).hexdigest(),
                      refit_method='OLS with intercept, Python standard library', python_version=platform.python_version())
    encoded = {path: table_bytes(path, rows) for path, rows in outputs.items()}
    encoded[base / 'calibration_refit.json'] = (json.dumps(provenance, indent=2) + '\n').encode()
    return encoded, summary


def run(families):
    prepared = [prepare(family, *FAMILIES[family]) for family in families]
    for outputs, summary in prepared:
        for path, content in outputs.items():
            path.write_bytes(content)
        print(f"{summary['series'].upper()}: T={summary['temperature_K']} K, "
              f"MAE={summary['calibrated_test_MAE_mean']:.6f}, "
              f"improved={summary['n_improved_vs_raw']}/100, "
              f"mean MAE reduction={summary['MAE_reduction_percent']:.3f}%")


if __name__ == '__main__':
    run(FAMILIES)
