"""Independent count-score/gradient checks and generating-truth leakage replay."""
import argparse
import hashlib
import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.provenance import start_run, finish_run


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


scalar = module('scalar', 'depth_reference_estimation_audit.py')
method = module('method', 'reference_depth_softmax.py')


def independent_gradient(theta, Y, t, scale, penalty):
    J = Y.shape[1]
    b = list(theta[:J-1])+[0j]
    q = list(theta[J-1:])+[0j]
    centre = sum(q)/J
    gb, gq = [0j]*J, [0j]*J
    total = float(Y.sum())
    for i, row in enumerate(Y):
        weights = [np.exp(z+t[i]*s) for z, s in zip(b, q)]
        normalizer = sum(weights)
        L = float(sum(row))
        for j, y in enumerate(row):
            residual = (L*weights[j]/normalizer-float(y))/total
            gb[j] += residual
            gq[j] += t[i]*residual
    gq = [v+penalty*(s-centre)/scale**2 for v, s in zip(gq, q)]
    return np.array(gb[:-1]+gq[:-1])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    panel = a.task/'panel_v1'
    source = a.task.parent/'depth_reference_estimation_20261109/panel_v2'
    run = start_run(a.out, dict(script=Path(__file__), design=a.task/'DESIGN.md',
        model=Path(method.__file__), scalar=Path(scalar.__file__), choices=panel/'choices.json',
        frozen=panel/'PREDICTIONS_FROZEN.json', data=source/'observations_and_targets.npz'), {},
        dict(scope='Numerical and leakage audit, no biological validation'), mode='reference_depth_softmax_audit')
    try:
        data = np.load(source/'observations_and_targets.npz')
        frozen = json.loads((panel/'PREDICTIONS_FROZEN.json').read_text())
        assert all(hashlib.sha256((panel/name).read_bytes()).hexdigest() == h for name, h in frozen.items())
        choices = json.loads((panel/'choices.json').read_text())['choices']
        frame = pd.read_csv(panel/'contrasts.tsv', sep='\t', keep_default_na=False)
        assert len(choices) == 48 and len(frame) == 480
        numerical, refinements, changed_targets = [], [], {}
        for choice in choices:
            condition, rep = choice['condition'], choice['rep']
            key = f'{condition}__{rep}'
            saved = np.load(panel/(key+'_predictions.npz'))
            Y = data[key+'__Y'].astype(float)
            H = data[key+'__H']
            L = Y.sum(axis=1)
            query = data[key+'__query_exposure']
            origin = float(saved['origin'])
            x = np.log(L)-origin
            scale = float(np.sqrt(np.mean(x*x)))
            budget = 256*Y.size*np.finfo(float).eps
            labels = ['pooled']+[f'ridge_{k}' for k in range(4, -5, -1)]+['unpenalized']
            scores, score_refinement = [], []
            for order, label in enumerate(labels):
                b, gamma = saved[label+'__b'], saved[label+'__gamma']
                reference = scalar.scalar_rates(b, gamma, L, origin)
                querymean = scalar.scalar_rates(b, gamma, query, origin).sum(axis=0)
                holdout = scalar.scalar_rates(b, gamma, H.sum(axis=1), origin)
                assert np.max(np.abs(querymean/saved[label+'__prediction']-1)) <= budget
                value = 2*math.fsum(float(m)-int(y)+(int(y)*math.log(int(y)/float(m)) if y else 0.)
                    for y, m in zip(H.flat, holdout.flat))/int(H.sum())
                assert abs(value-float(saved[label+'__score'])) <= budget
                scores.append((float(saved[label+'__score']), order, label))
                # Independent score equations in physical log-depth units.
                if label != 'pooled':
                    record = next(v for v in choice['solver_records'] if v['method'] == label)
                    penalty = record['fine']['penalty']
                    gb = [math.fsum(float(reference[i, j]-Y[i, j]) for i in range(len(Y)))/Y.sum() for j in range(40)]
                    gg = [math.fsum(float(x[i]*(reference[i, j]-Y[i, j])) for i in range(len(Y)))/Y.sum()
                          +penalty*(float(gamma[j])-math.fsum(gamma)/40) for j in range(40)]
                    assert max(abs(v) for v in gb[:-1]+[v/scale for v in gg[:-1]]) <= budget
                    assert record['fine']['min_Hessian_eigenvalue'] > 0
                    # Refit coarse controls to resolve whether selection overlaps numerical refinement.
                    cb, cg, _, _ = method.fit(Y, penalty, 1e-8)
                    coarse_score = method.score(H, method.means(cb, cg, H.sum(axis=1).astype(float), origin))
                    score_refinement.append(abs(coarse_score-float(saved[label+'__score'])))
                    if condition == 'supported' and rep == 0 and label == 'ridge_0':
                        theta = np.r_[b[:-1], gamma[:-1]*scale]
                        actual = method.objective(theta, Y, x/scale, scale, penalty, True)
                        assert np.max(np.abs(independent_gradient(theta.astype(complex), Y, x/scale, scale, penalty).real-actual[1])) <= budget
                        directions = [np.sin(np.arange(78)+1), np.cos(np.arange(78)+1), np.ones(78)]
                        errors = []
                        for direction in directions:
                            for step in (1e-15, 1e-25):
                                derivative = independent_gradient(theta+step*1j*direction, Y, x/scale, scale, penalty).imag/step
                                error = float(np.max(np.abs(derivative-actual[2]@direction)))
                                assert error <= budget
                                errors.append(error)
                        numerical.append(dict(independent_Hessian_directional_errors=errors, budget=budget))
                else:
                    score_refinement.append(0.)
            assert min(scores)[2] == choice['selected']
            ordered = sorted(scores)
            gap = ordered[1][0]-ordered[0][0]
            bound = score_refinement[ordered[0][1]]+score_refinement[ordered[1][1]]
            refinements.append(dict(condition=condition, rep=rep, selected=choice['selected'], score_gap=gap,
                coarse_fine_score_discrepancy=bound, selection_overlaps_refinement=bool(gap <= bound)))
            previous = np.load(source/(key+'_predictions.npz'))
            expected = dict(pooled=saved['pooled__prediction'], unpenalized=saved['unpenalized__prediction'],
                selected=saved[choice['selected']+'__prediction'], quintile=previous['estimated'], current=previous['current'])
            for event in ('null', 'gain'):
                target = data[key+'__'+event]
                effect = [.6*math.log(1.5) if event == 'gain' and j < 10 else 0. for j in range(40)]
                for label, E in expected.items():
                    values = [math.log(float(t)/float(e)) for t, e in zip(target, E)]
                    ordered_values = sorted(values)
                    centre = .5*(ordered_values[19]+ordered_values[20])
                    error = [v-centre-e for v, e in zip(values, effect)]
                    mse = math.fsum(v*v for v in error)/40
                    region = math.fsum(v-centre for v in values[:10])/10
                    selected = frame[(frame.condition == condition) & (frame.rep == rep) & (frame.event == event) & (frame.method == label)]
                    assert len(selected) == 1
                    row = selected.iloc[0]
                    assert abs(row.mse-mse) <= budget*max(1., mse)
                    assert abs(row.first_region-region) <= budget*max(1., abs(region))
                    assert abs(row.max_abs_error-max(abs(v) for v in error)) <= budget*max(1., max(abs(v) for v in error))
                changed_targets[key+'__'+event] = target*np.exp(.4*np.sin(np.arange(40)))
        # Anchor choice must not change an identified ridge prediction.
        Y = data['supported__0__Y'].astype(float)
        perm = np.arange(40)[::-1]
        selected = choices[0]
        penalty = next(v['fine']['penalty'] for v in selected['solver_records'] if v['method'] == 'ridge_0')
        b, g, origin, _ = method.fit(Y, penalty, 1e-10)
        bp, gp, op, _ = method.fit(Y[:, perm], penalty, 1e-10)
        query = data['supported__0__query_exposure']
        prediction = method.means(b, g, query, origin)
        permuted = method.means(bp, gp, query, op)[:, perm]
        permutation_error = float(np.max(np.abs(np.log(prediction/permuted))))
        assert permutation_error <= 256*Y.size*np.finfo(float).eps
        perturb = a.out/'perturbed_truth.npz'
        np.savez_compressed(perturb, **changed_targets)
        replay = a.out/'truth_replay'
        with (a.out/'replay_output.txt').open('w') as handle:
            subprocess.run([sys.executable, '-I', str(Path(method.__file__)), '--out', str(replay),
                '--task', str(a.task), '--truth', str(perturb)], stdout=handle, stderr=subprocess.STDOUT, check=True)
        for path in panel.glob('*_predictions.npz'):
            first, second = np.load(path), np.load(replay/path.name)
            assert first.files == second.files
            assert all(np.array_equal(first[name], second[name]) for name in first.files)
        for name in ('holdout_scores.tsv', 'choices.json'):
            assert (panel/name).read_bytes() == (replay/name).read_bytes()
        assert (panel/'contrasts.tsv').read_bytes() != (replay/'contrasts.tsv').read_bytes()
        pd.DataFrame(refinements).to_csv(a.out/'selection_refinement.tsv', sep='\t', index=False)
        counts = pd.DataFrame(refinements).groupby(['condition', 'selected']).size().rename('cases')
        counts.to_csv(a.out/'selected_models.tsv', sep='\t')
        (a.out/'report.json').write_text(json.dumps(dict(cases=48, scalar_prediction_bins=21120,
            scalar_heldout_scores=528, scalar_target_metrics=1440, independent_Hessian=numerical,
            anchor_permutation_error=permutation_error,
            selections_overlapping_refinement=int(sum(v['selection_overlaps_refinement'] for v in refinements)),
            prediction_truth_replay_exact=True, fitted_choices_truth_replay_exact=True,
            biological_specimens=0, production_changed=False, scientific_criteria_changed=False), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, run, error=error)
        raise
    finish_run(a.out, run, backend='CPU independent scalar math, score equations, Hessian and complete truth replay')


if __name__ == '__main__':
    main()
