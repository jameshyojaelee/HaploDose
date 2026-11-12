"""Check whether the zero-SD stationary point hides a better continuous null."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import expit
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    path = Path(__file__).with_name('nb_state_guard.py')
    spec = importlib.util.spec_from_file_location('boundary_nb_state', path)
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    inputs = dict(script=Path(__file__), model=path, design=args.task/'BOUNDARY_DESIGN.md')
    for seed in (20261107, 20261109):
        inputs[f'data{seed}'] = args.task/f'seed_{seed}_v1/proposal_observations.npz'
        inputs[f'fit{seed}'] = args.task/f'seed_{seed}_v1/report.json'
    rec = start_run(args.out, inputs, {}, dict(scope='Numerical optimum check, not scientific validation'),
        mode='gene_state_null_boundary')
    try:
        records = []
        for seed in (20261107, 20261109):
            with np.load(inputs[f'data{seed}']) as z:a = {k:z[k] for k in z.files}
            old = json.loads(inputs[f'fit{seed}'].read_text())
            for row in old['comparisons']:
                O, E, r, bins = a['counts'], a['expected'], a['size'], a['bins']
                if row['arm'] == 'moment_bin':
                    rr = np.stack([(E[:, bins==j]**2/r[:, bins==j]).sum(axis=1) for j in np.unique(bins)], axis=1)
                    O = np.stack([O[:, bins==j].sum(axis=1) for j in np.unique(bins)], axis=1)
                    E = np.stack([E[:, bins==j].sum(axis=1) for j in np.unique(bins)], axis=1)
                    r = E**2/rr
                fit = row['quadrature'][-1]
                mean = fit['null'][0]
                p = expit(np.log(E)+mean-np.log(r))
                score = (O-(O+r)*p).sum(axis=1)
                hessian = -((O+r)*p*(1-p)).sum(axis=1)
                # Derivative of negative mean log marginal likelihood in v=SD^2.
                derivative = -.5*float(np.mean(score**2+hessian))
                spread = float(np.std(np.log((O.sum(axis=1)+.5)/(E.sum(axis=1)+.5))))
                grid = [0.]+[spread*2.**k for k in range(-12, 4)]
                studies = []
                for sd in grid:
                    values = [model.objective(np.array([mean, sd]), O, E, r, order, False)[0]
                        for order in (24, 48, 96)]
                    studies.append(dict(sd=sd, objectives=values))
                best = min(studies, key=lambda x:x['objectives'][-1])
                operations = 1024*O.size*96+256
                gamma = operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)
                scale = max(1., abs(fit['null_objective']), abs(best['objectives'][-1]))
                discrepancy = abs(best['objectives'][-1]-best['objectives'][-2])
                improvement = fit['null_objective']-best['objectives'][-1]
                resolved = improvement > 2*discrepancy+gamma*scale
                records.append(dict(seed=seed, arm=row['arm'], null_variance_derivative_at_zero=derivative,
                    prior_null_objective=fit['null_objective'], prior_BIC=fit['BIC_advantage'],
                    grid=studies, best_fixed_mean=best, objective_improvement=improvement,
                    quadrature_discrepancy=discrepancy, rounding_budget=gamma*scale,
                    missed_better_null=bool(resolved),
                    BIC_upper_using_better_null=float(fit['BIC_advantage']-2*len(O)*improvement)))
        write_json(dict(records=records, confirmation_or_truth_used=False,
            frozen_results_changed=False), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU independently parameterized variance-boundary and fixed-mean grid')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':main()
