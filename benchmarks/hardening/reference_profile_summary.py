"""Descriptive frozen profile regression results; no prospective inference."""
import argparse
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    source = Path(__file__).with_name('performance_candidates.py')
    inputs = dict(script=Path(__file__), partition=source,
        scores=args.task/'regression_v1/scores.tsv', design=args.task/'REGRESSION_DESIGN.md',
        locked_scores=args.task/'prediction_v3/scores.tsv',
        identities=args.task/'regression_v1/predictions_locked.json')
    rec = start_run(args.out, inputs, dict(methods=['pseudobulk', 'pseudobulk_sm5', 'pseudobulk_mdl']),
        dict(scope='Previously exposed regression-only cell lines; shared references, no population CI or superiority test'),
        mode='reference_profile_summary')
    try:
        spec = importlib.util.spec_from_file_location('frozen_partition', source)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        checks = []
        for values in ([0., .1, -.1, 2., 2.1, 1.9], [0., 1., -1., .5, -.5, 1.]):
            y = np.array(values)
            v = np.array([.1, .2, .3, .15, .25, .2])
            fitted, parts = module.partition(y, v, np.r_[True, np.zeros(len(y)-1, bool)])
            penalty = 3*np.log(len(y))
            objective = float(np.sum((y-fitted)**2/v)+penalty*(len(parts)-1))
            candidates = []
            for cuts in itertools.product((False, True), repeat=len(y)-1):
                indices = np.r_[0, np.flatnonzero(cuts)+1, len(y)]
                value = penalty*(len(indices)-2)
                for lo, hi in zip(indices[:-1], indices[1:]):
                    mean = np.average(y[lo:hi], weights=1/v[lo:hi])
                    value += float(np.sum((y[lo:hi]-mean)**2/v[lo:hi]))
                candidates.append(value)
            independent = min(candidates)
            operations = 256*len(y)**2
            eps = np.finfo(float).eps
            budget = operations*eps/(1-operations*eps)*(1+max(map(abs, candidates)))
            assert abs(objective-independent) <= budget
            checks.append(dict(objective=objective, exhaustive_objective=independent,
                               rounding_budget=budget, partitions=len(candidates)))
        scores = pd.read_csv(inputs['scores'], sep='\t', float_precision='round_trip')
        old = pd.read_csv(inputs['locked_scores'], sep='\t', float_precision='round_trip')
        assert not scores.duplicated(['case', 'method']).any()
        matrix = scores.pivot(index='case', columns='method', values='pearson')
        assert matrix.shape == (181, 3) and np.isfinite(matrix.to_numpy()).all()
        assert scores.groupby('case').common_genes.nunique().eq(1).all()
        assert scores.groupby('case').windows.nunique().eq(1).all()
        locked_comparisons = 0
        for row in old[old.method.isin(matrix.columns)].itertuples():
            assert row.pearson == matrix.loc[row.case, row.method]
            locked_comparisons += 1
        rows = []
        for method in ('pseudobulk', 'pseudobulk_sm5'):
            delta = matrix.pseudobulk_mdl-matrix[method]
            rows.append(dict(comparator=method, cell_lines=len(delta), wins=int((delta > 0).sum()),
                ties=int((delta == 0).sum()), losses=int((delta < 0).sum()),
                mean_paired_change=float(delta.mean()), median_paired_change=float(delta.median()),
                worst_case=str(delta.idxmin()), worst_change=float(delta.min()),
                best_case=str(delta.idxmax()), best_change=float(delta.max())))
            matrix['change_vs_'+method] = delta
        identities = json.loads(inputs['identities'].read_text())['cases']
        assert len(identities) == 181 and len({r['case'] for r in identities}) == 181
        matrix.to_csv(args.out/'all_paired_profiles.tsv', sep='\t')
        write_json(dict(cell_lines=181, methods=scores.groupby('method').pearson.agg(['mean', 'median']).to_dict('index'),
            paired=rows, locked_criterion_values_identical=locked_comparisons,
            exhaustive_partition_checks=checks, shared_reference_dependence=True,
            population_superiority_established=False, independent_evaluation=False,
            production_changed=False), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU exhaustive small partition references and descriptive paired regression')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
