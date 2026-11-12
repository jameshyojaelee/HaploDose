"""Reference-only slope estimation and conditional RNA depth transport challenges."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.clones import reference_depth_slopes
from haplodose.provenance import start_run, finish_run

spec = importlib.util.spec_from_file_location('anchor', Path(__file__).with_name('depth_reference_anchor.py'))
anchor = importlib.util.module_from_spec(spec)
spec.loader.exec_module(anchor)


def predict(Y, gamma, query, tolerance):
    depth = Y.sum(axis=1).astype(float)
    mu = Y.sum(axis=0).astype(float) / Y.sum()
    origin = float(np.log(depth).mean())
    b, solver = anchor.calibrate(mu, gamma, depth, origin, tolerance)
    expected = query[:, None] * anchor.rates(b, gamma, query, origin)
    return expected, b, origin, solver


def generate(rng, n, limits, gamma):
    exposure = np.exp(rng.uniform(*np.log(limits), n))
    mean = exposure[:, None] * anchor.rates(np.zeros(40), gamma, exposure, math.log(9000.))
    Y = rng.poisson(mean * rng.gamma(20., .05, mean.shape)).astype(np.int64)
    if not (Y.sum(axis=1) > 0).all() or not (Y.sum(axis=0) > 0).all():
        raise ArithmeticError('Zero count support; no invented floor')
    return exposure, Y


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    run = start_run(a.out, dict(script=Path(__file__), design=a.task/'DESIGN.md',
        anchor=Path(anchor.__file__)), {}, dict(scope='Conditional raw RNA simulations; biological n=0'),
        mode='depth_reference_estimation')
    try:
        slope = np.r_[np.full(10, -.5), np.zeros(20), np.full(10, .5)]
        rows, checks, scores, data = [], [], [], {}
        for rep in range(8):
            for index, condition in enumerate(('supported', 'extrapolated', 'narrow', 'small', 'zero', 'confounded')):
                rng = np.random.default_rng(np.random.SeedSequence([2026110900+rep, index]))
                refgamma = np.zeros(40) if condition == 'zero' else slope.copy()
                querygamma = refgamma.copy()
                if condition == 'confounded':
                    refgamma[:10] += .5
                limits = (8500., 9500.) if condition == 'narrow' else (4500., 18000.)
                querylimits = limits if condition == 'supported' else (36000., 72000.)
                refexposure, Y = generate(rng, 30 if condition == 'small' else 200, limits, refgamma)
                holdexposure, H = generate(rng, 100, limits, refgamma)
                queryexposure = np.exp(rng.uniform(*np.log(querylimits), 100))
                groups = [np.arange(len(Y))]
                L = Y.sum(axis=1).astype(float)
                gamma = reference_depth_slopes(Y, L, groups, 30)
                gamma -= gamma.mean()
                mu = Y.sum(axis=0) / Y.sum()
                fine, b, origin, solver = predict(Y, gamma, queryexposure, 1e-12)
                coarse, _, _, coarse_solver = predict(Y, gamma, queryexposure, 1e-10)
                oracle, oracle_b, _, oracle_solver = predict(Y, refgamma, queryexposure, 1e-12)
                expected = dict(current=anchor.current_expected(mu, gamma, queryexposure,
                        float(np.log(queryexposure).mean())),
                    origin_only=anchor.current_expected(mu, gamma, queryexposure, origin),
                    estimated=fine.sum(axis=0), oracle=oracle.sum(axis=0))
                # Holdout observations are scored after fitting; they do not enter fitting.
                hold = H.sum(axis=1).astype(float)
                for method, g in [('pooled', np.zeros(40)), ('estimated', gamma), ('oracle', refgamma)]:
                    prediction, _, _, _ = predict(Y, g, hold, 1e-12)
                    positive = H > 0
                    term = prediction-H
                    term[positive] += H[positive]*np.log(H[positive]/prediction[positive])
                    scores.append(dict(rep=rep, condition=condition, method=method,
                        deviance_per_count=float(2*term.sum()/H.sum())))
                # Mutation and permutation checks use fit-only inputs; no regenerated truth.
                altered = H[::-1] * 7
                joint = np.vstack([Y, altered])
                g2 = reference_depth_slopes(joint, joint.sum(axis=1), groups, 30)
                g2 -= g2.mean()
                assert np.array_equal(gamma, g2)
                p2, _, _, _ = predict(joint[:len(Y)], g2, queryexposure, 1e-12)
                assert np.array_equal(p2, fine)
                perm = np.arange(len(Y))[::-1]
                gp = reference_depth_slopes(Y[perm], L[perm], groups, 30)
                gp -= gp.mean()
                eps = np.finfo(float).eps
                operations = 256*Y.size
                budget = operations*eps/(1-operations*eps)
                assert np.max(np.abs(gp-gamma)) <= budget
                assert np.max(np.abs(fine.sum(axis=1)-queryexposure)/queryexposure) <= budget
                key = f'{condition}__{rep}'
                # Freeze predictors and scores before constructing query target means.
                np.savez_compressed(a.out/(key+'_predictions.npz'), **expected, gamma=gamma,
                    b=b, oracle_b=oracle_b, origin=np.array(origin))
                for event in ('null', 'gain'):
                    effect = np.r_[np.full(10, .6*math.log(1.5)), np.zeros(30)] if event == 'gain' else np.zeros(40)
                    target = (queryexposure[:, None]*anchor.rates(effect, querygamma,
                        queryexposure, math.log(9000.))).sum(axis=0)
                    data[key+'__'+event] = target
                    for method, E in expected.items():
                        contrast = np.log(target/E)
                        contrast -= np.median(contrast)
                        error = contrast-effect
                        rows.append(dict(rep=rep, condition=condition, event=event, method=method,
                            mse=float(np.mean(error**2)), max_abs_error=float(np.max(np.abs(error))),
                            first_region=float(contrast[:10].mean())))
                checks.append(dict(rep=rep, condition=condition, solver=solver, coarse_solver=coarse_solver,
                    oracle_solver=oracle_solver, observed_reference_min=float(L.min()),
                    observed_reference_max=float(L.max()), query_min=float(queryexposure.min()),
                    query_max=float(queryexposure.max()), gamma_mse=float(np.mean((gamma-(refgamma-refgamma.mean()))**2)),
                    no_slope_support=bool(np.all(gamma == 0) and condition == 'small'),
                    max_permutation_change=float(np.max(np.abs(gp-gamma))), numerical_budget=budget,
                    max_refinement_change=float(np.max(np.abs(np.log(coarse/fine))))))
                for name, value in dict(Y=Y, H=H, reference_exposure=refexposure,
                        heldout_exposure=holdexposure, query_exposure=queryexposure,
                        reference_gamma=refgamma, query_gamma=querygamma).items():
                    data[key+'__'+name] = value
        pd.DataFrame(rows).to_csv(a.out/'contrasts.tsv', sep='\t', index=False)
        pd.DataFrame(scores).to_csv(a.out/'holdout_scores.tsv', sep='\t', index=False)
        np.savez_compressed(a.out/'observations_and_targets.npz', **data)
        (a.out/'report.json').write_text(json.dumps(dict(checks=checks, biological_specimens=0,
            simulation_replicates=8, production_changed=False, scientific_criteria_changed=False), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, run, error=error)
        raise
    finish_run(a.out, run, backend='CPU SciPy and installed reference_depth_slopes/expected')


if __name__ == '__main__':
    main()
