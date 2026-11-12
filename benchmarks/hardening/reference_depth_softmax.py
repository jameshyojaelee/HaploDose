"""Continuous fit-only RNA depth regression with heldout complexity selection."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize, root
from scipy.special import logsumexp
from haplodose.provenance import start_run, finish_run


def objective(theta, Y, t, scale, penalty, hessian=False):
    J = Y.shape[1]
    b = np.r_[theta[:J-1], 0.]
    q = np.r_[theta[J-1:], 0.]
    z = b[None]+t[:, None]*q[None]
    total = Y.sum()
    L = Y.sum(axis=1)
    mu = Y.sum(axis=0)/total
    baseline = np.log(mu)
    normalizer = logsumexp(z, axis=1)
    p = np.exp(z-normalizer[:, None])
    gamma = q/scale
    centred = gamma-gamma.mean()
    loss = (L@(normalizer-logsumexp(baseline))-np.sum(Y*(z-baseline)))/total
    loss += .5*penalty*(centred@centred)
    residual = (L[:, None]*p-Y)/total
    gb = residual.sum(axis=0)
    gq = (t[:, None]*residual).sum(axis=0)+penalty*centred/scale
    gradient = np.r_[gb[:-1], gq[:-1]]
    if not hessian:
        return float(loss), gradient
    blocks = []
    for power in range(3):
        w = L*t**power/total
        matrix = np.diag(w@p)-np.einsum('i,ij,ik->jk', w, p, p)
        if power == 2:
            matrix += penalty/scale**2*(np.eye(J)-np.ones((J, J))/J)
        blocks.append(matrix[:-1, :-1])
    H = np.block([[blocks[0], blocks[1]], [blocks[1], blocks[2]]])
    return float(loss), gradient, H


def scalar_objective(theta, Y, t, scale, penalty):
    J = Y.shape[1]
    b = np.r_[theta[:J-1], 0.]
    q = np.r_[theta[J-1:], 0.]
    total = 0j
    for i, row in enumerate(Y):
        logits = b+t[i]*q
        normalizer = np.log(sum(np.exp(v) for v in logits))
        total += sum(float(y)*(normalizer-v) for y, v in zip(row, logits))
    gamma = q/scale
    centred = gamma-sum(gamma)/J
    return total/Y.sum()+.5*penalty*sum(v*v for v in centred)


def fit(Y, penalty, tolerance):
    L = Y.sum(axis=1).astype(float)
    origin = float(np.log(L).mean())
    x = np.log(L)-origin
    scale = float(np.sqrt(np.mean(x*x)))
    if scale <= 0:
        raise ValueError('No identified reference depth variation')
    t = x/scale
    mu = Y.sum(axis=0)/Y.sum()
    start = np.r_[np.log(mu[:-1]/mu[-1]), np.zeros(Y.shape[1]-1)]
    optimization = minimize(objective, start, args=(Y, t, scale, penalty), jac=True,
        method='L-BFGS-B', options=dict(gtol=tolerance, ftol=0., maxiter=4000, maxls=50))
    polished = root(lambda z: objective(z, Y, t, scale, penalty)[1], optimization.x,
        jac=lambda z: objective(z, Y, t, scale, penalty, True)[2], options=dict(xtol=tolerance))
    value, gradient, H = objective(polished.x, Y, t, scale, penalty, True)
    operations = 128*Y.size
    unit = np.finfo(float).eps
    budget = operations*unit/(1-operations*unit)
    norm = float(np.max(np.abs(gradient)))
    if norm > max(tolerance, budget):
        raise ArithmeticError(f'Unresolved regression stationarity: {norm} > {max(tolerance, budget)}')
    J = Y.shape[1]
    b = np.r_[polished.x[:J-1], 0.]
    gamma = np.r_[polished.x[J-1:], 0.]/scale
    return b, gamma, origin, dict(penalty=float(penalty), tolerance=tolerance,
        optimizer_success=bool(optimization.success), optimizer_message=str(optimization.message),
        root_success=bool(polished.success), root_message=str(polished.message),
        gradient_max=norm, numerical_budget=budget, objective=value,
        min_Hessian_eigenvalue=float(np.linalg.eigvalsh(H).min()))


def means(b, gamma, L, origin):
    z = b[None]+(np.log(L)-origin)[:, None]*gamma[None]
    return L[:, None]*np.exp(z-logsumexp(z, axis=1, keepdims=True))


def score(Y, prediction):
    positive = Y > 0
    term = prediction-Y
    term[positive] += Y[positive]*np.log(Y[positive]/prediction[positive])
    return float(2*term.sum()/Y.sum())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--truth', type=Path, help='Evaluation-only source for leakage replay')
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    panel = a.task.parent/'depth_reference_estimation_20261109/panel_v2'
    truth_path = a.truth or panel/'observations_and_targets.npz'
    run = start_run(a.out, dict(script=Path(__file__), design=a.task/'DESIGN.md',
        observations=panel/'observations_and_targets.npz', truth=truth_path), {},
        dict(scope='Development RNA regression, biological n=0'), mode='reference_depth_softmax')
    try:
        data = np.load(panel/'observations_and_targets.npz')
        choices, checks, all_scores = [], [], []
        for condition in ('supported', 'extrapolated', 'narrow', 'small', 'zero', 'confounded'):
            for rep in range(8):
                key = f'{condition}__{rep}'
                Y = data[key+'__Y'].astype(float)
                H = data[key+'__H']
                query = data[key+'__query_exposure']
                L = Y.sum(axis=1)
                mu = Y.sum(axis=0)/Y.sum()
                origin = float(np.log(L).mean())
                x = np.log(L)-origin
                curvature = float(np.sum(L*x*x)/L.sum()*(1-mu@mu)/len(mu))
                candidates = [('pooled', None)]+[(f'ridge_{k}', curvature*10.**k) for k in range(4, -5, -1)]+[('unpenalized', 0.)]
                stored, records, rankings = {}, [], []
                for order, (label, penalty) in enumerate(candidates):
                    if penalty is None:
                        b, gamma = np.log(mu), np.zeros(40)
                    else:
                        coarse_b, coarse_g, _, coarse_record = fit(Y, penalty, 1e-8)
                        b, gamma, _, record = fit(Y, penalty, 1e-10)
                        fine_mean = means(b, gamma, query, origin)
                        coarse_mean = means(coarse_b, coarse_g, query, origin)
                        records.append(dict(method=label, coarse=coarse_record, fine=record,
                            max_refinement_log_prediction=float(np.max(np.abs(np.log(coarse_mean/fine_mean))))))
                        if rep == 0 and condition == 'supported' and label == 'ridge_0':
                            scale = float(np.sqrt(np.mean(x*x)))
                            theta = np.r_[np.log(mu[:-1]/mu[-1]), np.linspace(-.2, .2, 39)]
                            analytic = objective(theta, Y, x/scale, scale, penalty)[1]
                            discrepancy = []
                            for step in (1e-15, 1e-25):
                                numeric = []
                                for j in range(len(theta)):
                                    point = theta.astype(complex)
                                    point[j] += step*1j
                                    numeric.append(scalar_objective(point, Y, x/scale, scale, penalty).imag/step)
                                discrepancy.append(float(np.max(np.abs(analytic-numeric))))
                            assert max(discrepancy) <= record['numerical_budget']
                            checks.append(dict(complex_gradient_discrepancy=discrepancy, budget=record['numerical_budget']))
                    predicted = means(b, gamma, query, origin).sum(axis=0)
                    heldout = means(b, gamma, H.sum(axis=1).astype(float), origin)
                    deviance = score(H, heldout)
                    ref = means(b, gamma, L, origin).sum(axis=0)/L.sum()
                    budget = 128*Y.size*np.finfo(float).eps
                    assert np.max(np.abs(ref-mu)) <= max(1e-10, budget)
                    assert np.max(np.abs(means(b+gamma*.7, gamma, query, origin+.7).sum(axis=0)/predicted-1)) <= budget
                    stored[label+'__b'] = b
                    stored[label+'__gamma'] = gamma
                    stored[label+'__prediction'] = predicted
                    stored[label+'__score'] = np.array(deviance)
                    rankings.append((deviance, order, label))
                    all_scores.append(dict(condition=condition, rep=rep, method=label, score=deviance))
                selected = min(rankings)[2]
                choices.append(dict(condition=condition, rep=rep, selected=selected, curvature=curvature,
                    solver_records=records))
                stored['origin'] = np.array(origin)
                np.savez_compressed(a.out/(key+'_predictions.npz'), **stored)
        pd.DataFrame(all_scores).to_csv(a.out/'holdout_scores.tsv', sep='\t', index=False)
        (a.out/'choices.json').write_text(json.dumps(dict(choices=choices, numerical_checks=checks), indent=2)+'\n')
        paths = sorted(a.out.glob('*_predictions.npz'))+[a.out/'holdout_scores.tsv', a.out/'choices.json']
        frozen = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
        (a.out/'PREDICTIONS_FROZEN.json').write_text(json.dumps(frozen, indent=2)+'\n')
        truth = np.load(truth_path)
        rows = []
        for choice in choices:
            condition, rep = choice['condition'], choice['rep']
            key = f'{condition}__{rep}'
            saved = np.load(a.out/(key+'_predictions.npz'))
            previous = np.load(panel/(key+'_predictions.npz'))
            expected = dict(pooled=saved['pooled__prediction'], unpenalized=saved['unpenalized__prediction'],
                selected=saved[choice['selected']+'__prediction'], quintile=previous['estimated'], current=previous['current'])
            for event in ('null', 'gain'):
                target = truth[key+'__'+event]
                effect = np.r_[np.full(10, .6*math.log(1.5)), np.zeros(30)] if event == 'gain' else np.zeros(40)
                for method, E in expected.items():
                    contrast = np.log(target/E)
                    contrast -= np.median(contrast)
                    error = contrast-effect
                    rows.append(dict(condition=condition, rep=rep, event=event, method=method,
                        mse=float(np.mean(error**2)), max_abs_error=float(np.max(np.abs(error))),
                        first_region=float(contrast[:10].mean())))
        pd.DataFrame(rows).to_csv(a.out/'contrasts.tsv', sep='\t', index=False)
        summary = pd.DataFrame(rows).groupby(['condition', 'event', 'method'], sort=False).agg(
            mean_mse=('mse', 'mean'), mean_first_region=('first_region', 'mean'), worst_abs_error=('max_abs_error', 'max'))
        summary.to_csv(a.out/'summary.tsv', sep='\t')
        assert all(hashlib.sha256((a.out/name).read_bytes()).hexdigest() == digest for name, digest in frozen.items())
        (a.out/'report.json').write_text(json.dumps(dict(cases=48, candidate_models=528,
            biological_specimens=0, production_changed=False, scientific_criteria_changed=False,
            target_scope='Expected normalized RNA, not DNA copies or clone identity'), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, run, error=error)
        raise
    finish_run(a.out, run, backend='CPU SciPy continuous softmax mean regression')


if __name__ == '__main__':
    main()
