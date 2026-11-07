"""Independent arithmetic recount of frozen recovery evidence; compute nodes only.

No scientific decision threshold is changed. Bound discrepancies are exported in
absolute units and float64 ULPs for review, rather than hidden by a new tolerance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, special, stats


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            h.update(block)
    return dict(path=str(path.resolve()), sha256=h.hexdigest(), bytes=path.stat().st_size)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def bound(k, n, alpha, direction):
    """Invert binomial probabilities directly; no beta quantile calculation."""
    require(0 <= k <= n and n > 0 and 0 < alpha < 1, 'Invalid binomial endpoint')
    if direction == 'lower':
        return 0. if k == 0 else float(optimize.brentq(
            lambda p: stats.binom.sf(k - 1, n, p) - alpha, 0., 1.,
            xtol=np.nextafter(0., 1.), rtol=4 * np.finfo(float).eps))
    require(direction == 'upper', 'Unknown bound direction')
    return 1. if k == n else float(optimize.brentq(
        lambda p: stats.binom.cdf(k, n, p) - alpha, 0., 1.,
        xtol=np.nextafter(0., 1.), rtol=4 * np.finfo(float).eps))


def discrepancy(saved, reference):
    require(np.isfinite([saved, reference]).all(), 'Missing/nonfinite saved bound')
    delta = float(saved - reference)
    return dict(saved=float(saved), reference=float(reference), difference=delta,
                absolute_difference=abs(delta),
                ulp_difference=abs(delta) / max(abs(float(np.spacing(reference))), np.nextafter(0., 1.)))


def jsonlines(path):
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    require(bool(rows), f'Empty required input: {path}')
    return rows


def clone_check(root):
    summary = json.loads((root / 'summary.json').read_text())
    design = summary['design']
    require(json.loads((root / 'manifest.json').read_text())['design'] == design,
            'Clone executed manifest and summary designs differ')
    require(design['alpha'] == .05 and design['evaluation_replicates_per_condition'] == 200
            and design['calibration_replicates_per_null'] == 199, 'Clone operating point changed')
    conditions = ['matched_nb', 'correlated_state', 'two_clones']
    methods = ['raw', 'NB_calibrated_gate', 'state_calibrated_gate']
    require(design['evaluation_conditions'] == conditions and design['methods'] == methods, 'Clone endpoint family changed')
    calibration, evaluation = jsonlines(root / 'calibration.jsonl'), jsonlines(root / 'evaluation.jsonl')
    require(len(calibration) == 398 and len(evaluation) == 600, 'Incomplete clone outputs')
    banks = {}
    for condition in conditions[:2]:
        rows = [r for r in calibration if r['condition'] == condition]
        require(len(rows) == 199 and {r['replicate'] for r in rows} == set(range(199)), 'Incomplete/duplicate calibration bank')
        values = np.array([r['score'] for r in rows], float)
        require(np.isfinite(values).all(), 'Nonfinite calibration score')
        banks[condition] = values
    rank_checks, endpoints = [], []
    require(len(summary['summary']) == 9, 'Missing clone summary endpoints')
    saved = {(r['condition'], r['method']): r for r in summary['summary']}
    require(len(saved) == 9, 'Duplicate clone summary endpoints')
    for condition in conditions:
        rows = [r for r in evaluation if r['condition'] == condition]
        require(len(rows) == 200 and {r['replicate'] for r in rows} == set(range(200)), 'Incomplete/duplicate evaluation condition')
        events = dict.fromkeys(methods, 0)
        for row in rows:
            require(np.isfinite(row['score']), 'Nonfinite evaluation score')
            expected = {'raw': row['raw_candidate_count'] > 1}
            for method, bank in zip(methods[1:], conditions[:2]):
                p = (1 + np.count_nonzero(banks[bank] >= row['score'])) / 200
                actual = row['methods'][method]['rank_p']
                require(actual == p, 'Saved rank differs from bank recount')
                expected[method] = bool(row['raw_candidate_count'] > 1 and row['score'] > 0 and p <= .05)
                rank_checks.append(dict(condition=condition, replicate=row['replicate'], method=method,
                                        rank_p=p, saved_rank_p=actual))
            for method in methods:
                require(row['methods'][method]['accepted_multigroup'] == expected[method], 'Saved acceptance decision changed')
                events[method] += int(expected[method])
        for method in methods:
            old = saved[condition, method]
            direction = 'lower' if condition == 'two_clones' else 'upper'
            require(old['n'] == 200 and old['events'] == events[method]
                    and old['rate'] == events[method] / 200
                    and old['joint_family9_direction'] == direction, 'Clone endpoint recount differs')
            endpoints.append(dict(condition=condition, method=method, events=events[method], n=200,
                                  direction=direction, tail_alpha=.05 / 9,
                                  bound=discrepancy(old['joint_family9_bound'], bound(events[method], 200, .05 / 9, direction))))
    return dict(endpoints=endpoints, ranks_recounted=len(rank_checks), ranks=rank_checks,
                scope='nine bounds conditional on realized calibration banks; no biological replication')


def uncertainty_check(root):
    saved = json.loads((root / 'summary.json').read_text())
    design = saved['design']
    require(design['coverage_target'] == .95 and design['family_alpha'] == .05
            and design['replicates_per_condition'] == 1000, 'Original coverage criterion changed')
    conditions = ['matched_homogeneous', 'matched_heterogeneous', 'correlated_heterogeneous']
    selections = ['prespecified', 'largest_absolute_mean_among_HMM_runs']
    methods = ['legacy_reported_se', 'simultaneous']
    require(design['conditions'] == conditions and design['intervals'] == selections, 'Coverage family changed')
    frame = pd.read_csv(root / 'intervals.tsv', sep='\t')
    require(set(frame.phase) <= {'coverage', 'power'}, 'Unexpected experiment phase')
    coverage = frame[frame.phase.eq('coverage')].copy()
    require(len(coverage) == 12000, 'Incomplete/extra primary coverage rows')
    require(np.isfinite(coverage.true_mu_I).all() and coverage.true_mu_I.eq(0).all(), 'Changed population-mu null target')
    families = {}
    for method in methods:
        require(len(saved['coverage_families'][method]) == 6, 'Incomplete six-endpoint summary')
        old = {(r['condition'], r['interval']): r for r in saved['coverage_families'][method]}
        require(len(old) == 6, 'Duplicate uncertainty endpoints')
        rows, logps = [], []
        for condition in conditions:
            for selection in selections:
                x = coverage[coverage.method.eq(method) & coverage.condition.eq(condition) & coverage.interval.eq(selection)]
                require(len(x) == 1000 and set(x.replicate) == set(range(1000)), 'Incomplete/duplicate coverage endpoint')
                require(not x[['lower', 'upper', 'true_mu_I']].isna().any().any()
                        and (x.lower <= x.upper).all(), 'Invalid required interval')
                covered = (x.lower <= x.true_mu_I) & (x.true_mu_I <= x.upper)
                require(np.array_equal(covered, x.covered), 'Saved coverage flags differ from true_mu target')
                k, n = int(covered.sum()), len(x)
                previous = old[condition, selection]
                require(previous['coverage']['successes'] == k and previous['coverage']['n'] == n
                        and previous['coverage']['fraction'] == k / n, 'Coverage count differs')
                # Direct log-binomial summation retains tails that CDF underflows.
                support = np.arange(k + 1)
                logp = float(special.logsumexp(special.gammaln(n + 1) - special.gammaln(support + 1)
                    - special.gammaln(n - support + 1) + support * np.log(.95) + (n - support) * np.log(.05)))
                logps.append(logp)
                comparisons = {}
                for label, confidence in [('coverage', .95), ('simultaneous_cp', 1 - .05 / 6)]:
                    tail = (1 - confidence) / 2
                    comparisons[label] = {name: discrepancy(previous[label][key], bound(k, n, tail, name))
                                          for name, key in [('lower', 'cp_low'), ('upper', 'cp_high')]}
                rows.append(dict(condition=condition, interval=selection, successes=k, n=n,
                                 bounds=comparisons, log_p_undercoverage=logp))
        order = np.argsort(logps, kind='stable')
        adjusted = np.empty(6)
        adjusted[order] = np.minimum(0., np.maximum.accumulate(np.asarray(logps)[order] + np.log(np.arange(6, 0, -1))))
        for row, lp in zip(rows, adjusted):
            previous = old[row['condition'], row['interval']]
            require(previous['falsified_undercoverage'] == bool(lp < np.log(.05)), 'Original six-endpoint decision differs')
            row.update(log_p_holm=float(lp), saved_log_p_holm=previous['log_p_holm'],
                       log_p_holm_difference=float(previous['log_p_holm'] - lp),
                       falsified_undercoverage=bool(lp < np.log(.05)))
        families[method] = rows
    return dict(coverage_families=families, secondary_power_rows_excluded_from_primary=int(frame.phase.eq('power').sum()),
                target='same population systematic mu_I=0', no_new_scientific_test=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--clone-root', type=Path, required=True)
    parser.add_argument('--uncertainty-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    roots, out = [args.clone_root.resolve(), args.uncertainty_root.resolve()], args.out.resolve()
    require(not out.exists() and not any(r == out or r in out.parents for r in roots), 'Output must be new and outside inputs')
    out.mkdir(parents=True)
    report = dict(status='failed', scientific_thresholds_unchanged=True,
                  arithmetic_policy='Report reference discrepancies without changing scientific decision thresholds')
    files = [roots[0] / name for name in ['manifest.json', 'summary.json', 'calibration.jsonl', 'evaluation.jsonl']]
    files += [roots[1] / name for name in ['summary.json', 'intervals.tsv', 'source_manifest.json']]
    try:
        report['immutable_inputs_before'] = [digest(p) for p in files]
        report['clones'] = clone_check(roots[0])
        report['uncertainty'] = uncertainty_check(roots[1])
        report['immutable_inputs_after'] = [digest(p) for p in files]
        require(report['immutable_inputs_before'] == report['immutable_inputs_after'], 'Input changed during audit')
        report['checker_source'] = digest(Path(__file__))
        report['status'] = 'completed_reference_recount'
    except Exception as exc:
        report['error'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
