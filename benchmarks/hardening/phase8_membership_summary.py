"""Locate loss of an exposed simulated contrast after an exact frozen replay."""
import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.special import xlogy
from scipy.stats import rankdata
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def auc(values, truth):
    """Descriptive rank separation, independently checked by pair enumeration."""
    a, b = values[truth == 1], values[truth == 0]
    ranks = rankdata(values)
    value = (ranks[truth == 1].sum() - len(a)*(len(a)+1)/2)/(len(a)*len(b))
    direct = (np.sum(a[:, None] > b) + .5*np.sum(a[:, None] == b))/(len(a)*len(b))
    assert value == direct
    return float(value)


def membership(R, truth):
    assert R.shape[0] == len(truth) and np.isfinite(R).all() and np.all(R >= 0)
    labels = R.argmax(axis=1)
    table = np.array([[np.sum((truth == t) & (labels == k))
                       for k in range(1, R.shape[1])] for t in (0, 1)])
    i, j = linear_sum_assignment(-table)
    return dict(mass=R.sum(axis=0).tolist(), map_table=table.tolist(),
                normal_maps=int(np.sum(labels == 0)),
                all_cell_map_correct=int(table[i, j].sum()),
                mean_entropy_nats=float(-xlogy(R, R).sum(axis=1).mean()),
                responsibility_truth_auc=[auc(R[:, k], truth) for k in range(1, R.shape[1])])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--replicate', type=int, choices=(1, 2, 4), default=1)
    args = parser.parse_args()
    seed = 202611098000+args.replicate
    trace = args.task/'trace_v1'
    assert json.loads((trace/'run_manifest.json').read_text())['status'] == 'complete'
    report = json.loads((trace/'report.json').read_text())
    assert len(report['exact_output_parity']) == 2
    assert all(r['exact'] and r['arrays'] and r['tables'] for r in report['exact_output_parity'])
    case = trace/'replay'/f'seed_{seed}'
    candidate = case/'pooled_prevalence'
    records = json.loads((trace/'stages/records.json').read_text())
    guards = json.loads((candidate/'guards/records.json').read_text())
    inputs = dict(script=Path(__file__), design=args.task/'SUMMARY_DESIGN.md',
                  replay_report=trace/'report.json', replay_manifest=trace/'run_manifest.json',
                  records=trace/'stages/records.json', guards=candidate/'guards/records.json',
                  generation=case/'generation.npz', bins=candidate/'fit/bins.tsv')
    for n, row in enumerate(records):
        inputs[f'stage_{n}'] = Path(row['arrays'])
    for row in guards:
        if row.get('guard_passed'):
            folder = Path(row['folder'])
            inputs[f"guard_{row['number']}_discovery"] = folder/'discovery.npz'
            inputs[f"guard_{row['number']}_prediction"] = folder/'prediction_96.npz'
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, inputs, dict(seed=seed, replicate=args.replicate, refits=0),
                    dict(scope='Exposed conditional mechanism description; no independent error or power estimate'),
                    mode='phase8_membership_summary')
    try:
        with np.load(inputs['generation']) as z:
            query_rows, truth = z['query_rows'], z['truth']
        truth_index = pd.Series(truth, index=query_rows)
        bins = pd.read_csv(inputs['bins'], sep='\t')
        target = bins.chromosome.eq('chr12').to_numpy() & (bins.end.to_numpy() <= 10000000)
        # The synthetic event is gene-defined. Keep boundary bins separate;
        # this is a profile-location diagnostic, not resolution-matched DNA scoring.
        outside = ~bins.chromosome.eq('chr12').to_numpy()
        summaries, comparisons, guard_rows = [], [], []
        previous = None
        for n, row in enumerate(records):
            with np.load(row['arrays']) as z:
                selected_truth = truth_index.loc[z['rows']].to_numpy()
                assert len(selected_truth) == len(truth) == 600
                sides = {}
                for side in ('before', 'after'):
                    dose = z['doses_'+side]
                    assert dose.shape[1] == len(bins)
                    effect = np.median(dose[:, target], axis=1) - np.median(dose[:, outside], axis=1)
                    pairs = []
                    for a, b in combinations(range(1, len(dose)), 2):
                        contrast = dose[a]-dose[b]
                        pairs.append(dict(a=a, b=b, target_effect=float(effect[a]-effect[b]),
                            identical=bool(np.array_equal(dose[a], dose[b])),
                            max_absolute_contrast=float(np.max(np.abs(contrast))),
                            outside_rms=float(np.sqrt(np.mean(contrast[outside]**2)))))
                    sides[side] = dict(groups=len(dose)-1, target_effect=effect.tolist(), pairs=pairs)
                    if 'R_'+side in z:
                        sides[side]['membership'] = membership(z['R_'+side], selected_truth)
                change = {name: float(np.max(np.abs(z[name+'_after']-z[name+'_before'])))
                          for name in ('base', 'offset', 'phi', 'gamma', 'U')}
                item = dict(index=n, stage=row['stage'], log=row['log'],
                            before=sides['before'], after=sides['after'], nuisance_max_change=change)
                if row['stage'] == 'em':
                    LL = z['LL_after']
                    item['pair_loglik_truth_auc'] = [dict(a=a, b=b, auc=auc(LL[:, a]-LL[:, b], selected_truth))
                        for a, b in combinations(range(1, LL.shape[1]), 2)]
                if previous is not None:
                    item['previous_output_profiles_exact'] = bool(np.array_equal(previous, z['doses_before']))
                previous = z['doses_after'].copy()
                summaries.append(item)
                for guard in guards:
                    if not guard.get('guard_passed') or 'R_after' not in z:
                        continue
                    folder = Path(guard['folder'])
                    with np.load(folder/'discovery.npz') as d, np.load(folder/'prediction_96.npz') as p:
                        indices = pd.Index(z['rows']).get_indexer(d['global_rows'])
                        assert np.all(indices >= 0)
                        t = selected_truth[indices]
                        for key in ('discovery_membership', 'joint_membership'):
                            q = p[key]
                            assert q.shape == t.shape
                            # Compare two finite probability vectors without fitting or changing inference.
                            for k in range(1, z['R_after'].shape[1]):
                                r = z['R_after'][indices, k]
                                centered_q, centered_r = q-q.mean(), r-r.mean()
                                norm = np.linalg.norm(centered_q)*np.linalg.norm(centered_r)
                                comparisons.append(dict(stage=n, guard=guard['number'], kind=key,
                                    component=k, cells=len(t),
                                    correlation=float(centered_q@centered_r/norm) if norm > 0 else None))
                            if n == 0:
                                guard_rows.append(dict(guard=guard['number'], source=guard['source'],
                                    kind=key, score=guard['score'], cells=len(t), truth_auc=auc(q, t)))
        assert len(summaries) == len(records)
        if args.replicate in (1, 2):
            assert any(r['after']['groups'] == 1 for r in summaries)
        assert summaries[-1]['after']['groups'] == {1:1, 2:1, 4:3}[args.replicate]
        write_json(dict(stages=summaries, guards=guard_rows, guard_stage_correlations=comparisons,
                        refits=0, biological_specimens=0, calibration_claim=False,
                        exact_replay_gate_passed=True), args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU independent saved-stage NumPy/SciPy diagnostics')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
