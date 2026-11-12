"""Recount conditional bootstrap intervals without fitting or quantile helpers."""
import argparse
import json
import math
import statistics
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.provenance import start_run, finish_run

spec = importlib.util.spec_from_file_location('slope_audit', Path(__file__).with_name('depth_reference_estimation_audit.py'))
slope_audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(slope_audit)


def centred(values):
    logarithms = [math.log(float(v)) for v in values]
    median = statistics.median(logarithms)
    return [v-median for v in logarithms]


def exact_upper(k, n, alpha):
    if k == n:
        return 1.
    low, high = 0., 1.
    for _ in range(80):
        p = .5*(low+high)
        cdf = math.fsum(math.comb(n, j)*p**j*(1-p)**(n-j) for j in range(k+1))
        if cdf > alpha:
            low = p
        else:
            high = p
    return .5*(low+high)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    folder = a.task/'uncertainty_v1'
    panel = a.task/'panel_v2'
    run = start_run(a.out, dict(script=Path(__file__), design=a.task/'UNCERTAINTY_DESIGN.md',
        report=folder/'report.json', intervals=folder/'intervals.tsv',
        observations=panel/'observations_and_targets.npz', holdout_scores=panel/'holdout_scores.tsv',
        scalar_helper=Path(slope_audit.__file__)), {},
        dict(scope='Independent bootstrap quantile, interval and binomial recount'),
        mode='depth_reference_uncertainty_audit')
    try:
        data = np.load(panel/'observations_and_targets.npz')
        frame = pd.read_csv(folder/'intervals.tsv', sep='\t', keep_default_na=False)
        report = json.loads((folder/'report.json').read_text())
        scores = pd.read_csv(panel/'holdout_scores.tsv', sep='\t', keep_default_na=False)
        assert len(frame) == 96 and report['bootstrap_samples'] == 199 and report['quantile_order'] == 190
        records = []
        for summary in report['summaries']:
            condition = summary['condition']
            null_failures, null_regions, gains = 0, 0, 0
            for rep in range(8):
                key = f'{condition}__{rep}'
                original = np.load(panel/(key+'_predictions.npz'))
                H = data[key+'__H']
                hold_depth = H.sum(axis=1)
                Y = data[key+'__Y']
                pooled = Y.sum(axis=0)/Y.sum()
                score_budget = 256*H.size*np.finfo(float).eps
                for method, b, gamma in [('pooled', np.log(pooled), np.zeros(40)),
                        ('estimated', original['b'], original['gamma']),
                        ('oracle', original['oracle_b'], data[key+'__reference_gamma'])]:
                    predicted = slope_audit.scalar_rates(b, gamma, hold_depth, float(original['origin']))
                    deviance = 2*math.fsum(float(m)-int(y)+(int(y)*math.log(int(y)/float(m)) if y else 0.)
                        for y, m in zip(H.flat, predicted.flat))/int(H.sum())
                    score = scores[(scores.condition == condition) & (scores.rep == rep) & (scores.method == method)]
                    assert len(score) == 1
                    assert abs(float(score.iloc[0].deviance_per_count)-deviance) <= score_budget
                samples = np.load(folder/(key+'_bootstrap.npz'))
                base = centred(original['estimated'])
                budget = 64*40*199*np.finfo(float).eps
                assert np.max(np.abs(np.array(base)-samples['centre'])) <= budget
                assert samples['samples'].shape == (199, 40)
                indices = samples['resampled_indices']
                assert indices.shape == (199, len(data[key+'__Y']))
                assert ((indices >= 0) & (indices < indices.shape[1])).all()
                differences = [[float(v)-b for v, b in zip(row, base)] for row in samples['samples']]
                sup = sorted(max(abs(v) for v in row) for row in differences)
                region = sorted(abs(math.fsum(row[:10])/10) for row in differences)
                critical, fixed = sup[189], region[189]
                for event in ('null', 'gain'):
                    target = centred(data[key+'__'+event])
                    difference = [t-b for t, b in zip(target, base)]
                    value = math.fsum(difference[:10])/10
                    selected = frame[(frame.condition == condition) & (frame.rep == rep) & (frame.event == event)]
                    assert len(selected) == 1
                    row = selected.iloc[0]
                    for actual, expected in [(row.whole_profile_band, critical), (row.fixed_region_band, fixed),
                            (row.region, value), (row.region_lower, value-fixed), (row.region_upper, value+fixed)]:
                        assert abs(float(actual)-expected) <= budget*max(1., abs(expected))
                    # Require identical decisions; do not hide a boundary discrepancy with tolerance.
                    assert bool(row.positive_region) == (value-fixed > 0)
                    assert bool(row.excludes_zero) == (abs(value) > fixed)
                    failure = event == 'null' and max(abs(v) for v in difference) > critical
                    assert bool(row.null_profile_uncovered) == failure
                    if event == 'null':
                        null_failures += int(failure)
                        null_regions += int(abs(value) > fixed)
                    else:
                        gains += int(value-fixed > 0)
                records.append(dict(condition=condition, rep=rep, profile_band=critical,
                    region_band=fixed, arithmetic_budget=budget))
            assert null_failures == summary['whole_profile_null_failures']
            assert null_regions == summary['fixed_region_null_exclusions']
            assert gains == summary['fixed_region_gain_positive']
            for alpha, name in [(.05, 'null_failure_upper95'), (.05/6, 'null_failure_upper95_simultaneous6')]:
                bound = exact_upper(null_failures, 8, alpha)
                assert abs(bound-summary[name]) <= 128*np.finfo(float).eps
        (a.out/'report.json').write_text(json.dumps(dict(cases=len(records), bootstrap_quantiles=96,
            interval_decisions=288, exact_binomial_bounds=12, heldout_scalar_deviances=144, checks=records,
            all_checks_passed=True, biological_specimens=0,
            interpretation='Numerical interval recount, not biological coverage validation'), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, run, error=error)
        raise
    finish_run(a.out, run, backend='CPU scalar math/statistics and exact binomial polynomial inversion')


if __name__ == '__main__':
    main()
