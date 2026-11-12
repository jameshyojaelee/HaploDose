"""Independent scalar recount of the estimated-slope depth transport panel."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.provenance import start_run, finish_run


def scalar_rates(b, gamma, depth, origin):
    result = []
    for L in depth:
        logits = [float(z)+float(g)*(math.log(float(L))-origin) for z, g in zip(b, gamma)]
        maximum = max(logits)
        weights = [math.exp(z-maximum) for z in logits]
        denominator = math.fsum(weights)
        result.append([float(L)*w/denominator for w in weights])
    return np.array(result)


def scalar_slopes(Y):
    if len(Y) < 60:
        return np.zeros(Y.shape[1])
    L = [int(sum(row)) for row in Y]
    x = np.array([math.log(v) for v in L])
    x -= math.fsum(x)/len(x)
    # Explicit linear order-statistic interpolation, not np.quantile.
    ordered = sorted(x)
    edges = []
    for quantile in (.2, .4, .6, .8):
        index = (len(x)-1)*quantile
        low = math.floor(index)
        edges.append(ordered[low]+(index-low)*(ordered[math.ceil(index)]-ordered[low]))
    members = [[i for i, value in enumerate(x) if sum(value > edge for edge in edges) == group]
               for group in range(5)]
    totals = [sum(L[i] for i in idx) for idx in members]
    average = [math.fsum(L[i]*float(x[i]) for i in idx)/total for idx, total in zip(members, totals)]
    xc = [value-math.fsum(average)/5 for value in average]
    denominator = math.fsum(v*v for v in xc)
    gamma = []
    for j in range(Y.shape[1]):
        logshare = [math.log((sum(int(Y[i, j]) for i in idx)+.5)/total)
                    for idx, total in zip(members, totals)]
        centre = math.fsum(logshare)/5
        slope = math.fsum(v*(q-centre) for v, q in zip(xc, logshare))/denominator
        gamma.append(min(1.5, max(-1.5, slope)))
    centre = math.fsum(gamma)/len(gamma)
    return np.array([g-centre for g in gamma])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    panel = a.task/'panel_v2'
    a.out.mkdir(exist_ok=False)
    run = start_run(a.out, dict(script=Path(__file__), design=a.task/'DESIGN.md',
        data=panel/'observations_and_targets.npz', contrasts=panel/'contrasts.tsv',
        report=panel/'report.json'), {}, dict(scope='Independent numerical audit; biological n=0'),
        mode='depth_reference_estimation_audit')
    try:
        data = np.load(panel/'observations_and_targets.npz')
        frame = pd.read_csv(panel/'contrasts.tsv', sep='\t', keep_default_na=False)
        checks = json.loads((panel/'report.json').read_text())['checks']
        assert len(checks) == 48 and len(frame) == 384
        records, paired = [], []
        for check in checks:
            condition, rep = check['condition'], check['rep']
            key = f'{condition}__{rep}'
            prediction = np.load(panel/(key+'_predictions.npz'))
            Y = data[key+'__Y']
            query = data[key+'__query_exposure']
            budget = check['numerical_budget']
            g = scalar_slopes(Y)
            error = float(np.max(np.abs(g-prediction['gamma'])))
            assert error <= budget
            pooled = Y.sum(axis=0)/Y.sum()
            depth = Y.sum(axis=1)
            b = prediction['b']
            origin = float(prediction['origin'])
            expected = scalar_rates(b, g, query, origin).sum(axis=0)
            assert np.max(np.abs(expected/prediction['estimated']-1)) <= budget
            ref = scalar_rates(b, g, depth, origin).sum(axis=0)/depth.sum()
            assert np.max(np.abs(ref-pooled)) <= budget
            oracle = scalar_rates(prediction['oracle_b'], data[key+'__reference_gamma'], query, origin).sum(axis=0)
            assert np.max(np.abs(oracle/prediction['oracle']-1)) <= budget
            for event in ('null', 'gain'):
                delta = .6*math.log(1.5) if event == 'gain' else 0.
                effect = np.array([delta if j < 10 else 0. for j in range(40)])
                truth = scalar_rates(effect, data[key+'__query_gamma'], query, math.log(9000.)).sum(axis=0)
                assert np.max(np.abs(truth/data[key+'__'+event]-1)) <= budget
                group = frame[(frame.condition == condition) & (frame.rep == rep) & (frame.event == event)]
                assert set(group.method) == {'current', 'origin_only', 'estimated', 'oracle'}
                for row in group.itertuples():
                    contrast = [math.log(float(t)/float(e)) for t, e in zip(truth, prediction[row.method])]
                    sorted_contrast = sorted(contrast)
                    centre = .5*(sorted_contrast[19]+sorted_contrast[20])
                    residual = [q-centre-v for q, v in zip(contrast, effect)]
                    mse = math.fsum(v*v for v in residual)/40
                    maximum = max(abs(v) for v in residual)
                    region = math.fsum(contrast[j]-centre for j in range(10))/10
                    assert abs(mse-row.mse) <= budget*max(1., mse)
                    assert abs(maximum-row.max_abs_error) <= budget*max(1., maximum)
                    assert abs(region-row.first_region) <= budget*max(1., abs(region))
                values = group.set_index('method')
                paired.append(dict(rep=rep, condition=condition, event=event,
                    estimated_minus_current_mse=float(values.loc['estimated', 'mse']-values.loc['current', 'mse'])))
            records.append(dict(condition=condition, rep=rep, slope_max_error=error,
                budget=budget, refinement=check['max_refinement_change']))
        summary = frame.groupby(['condition', 'event', 'method'], sort=False).agg(
            mean_mse=('mse', 'mean'), min_mse=('mse', 'min'), max_mse=('mse', 'max'),
            mean_first_region=('first_region', 'mean'), worst_abs_error=('max_abs_error', 'max'))
        summary.to_csv(a.out/'summary.tsv', sep='\t')
        pd.DataFrame(paired).to_csv(a.out/'paired_mse.tsv', sep='\t', index=False)
        (a.out/'report.json').write_text(json.dumps(dict(checks=records, scalar_metrics=1152,
            biological_specimens=0, cases=48, all_checks_passed=True,
            interpretation='Numerical agreement only; no production/error/power/biological claim'), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, run, error=error)
        raise
    finish_run(a.out, run, backend='CPU independent scalar math/order statistics')


if __name__ == '__main__':
    main()
