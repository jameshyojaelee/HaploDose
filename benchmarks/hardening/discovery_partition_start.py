"""A discovery-only initializer probe for the unchanged shared-state likelihood."""
import argparse
from fractions import Fraction
import importlib.util
from itertools import product
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import minimize

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def load(name):
    spec = importlib.util.spec_from_file_location('partition_start_'+name,
        Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def partition(values):
    """Exact search over sorted two-mean cuts; float SSE, not a clone test."""
    values = np.sort(np.asarray(values, float))
    if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all():
        raise ValueError('At least two finite discovery amplitudes required')
    candidates = []
    for cut in range(1, len(values)):
        means = [math.fsum(map(float, v))/len(v) for v in (values[:cut], values[cut:])]
        cost = math.fsum(float(x-m)**2 for v, m in zip((values[:cut], values[cut:]), means)
                         for x in v)
        candidates.append((cost, cut, means))
    cost, cut, means = min(candidates, key=lambda r: (r[0], r[1]))
    return dict(cut=cut, cells=len(values), proportion=(len(values)-cut)/len(values),
                delta=means[1]-means[0], variance=cost/len(values),
                means=means, objective=cost, scientific_decision=False)


def numerical_checks():
    x = [-3, -2, -1, 1, 2, 3]
    reference = []
    for allocation in product((0, 1), repeat=len(x)):
        groups = [[Fraction(v) for v, a in zip(x, allocation) if a == k] for k in (0, 1)]
        if any(not g for g in groups):
            continue
        cost = sum(sum((v-sum(g)/len(g))**2 for v in g) for g in groups)
        reference.append(cost)
    assert min(reference) == Fraction(4)
    result = partition(x)
    assert result['cut'] == 3 and result['objective'] == 4. and result['delta'] == 4.
    assert partition(x[::-1]) == result
    assert partition([v+1024 for v in x])['objective'] == 4.
    scaled = partition([v*8 for v in x])
    assert scaled['objective'] == 256. and scaled['delta'] == 32.
    p, delta, null_variance, variance = .25, math.log(2.), .5, .125
    mean = .5*(null_variance-variance)-math.log1p(p*math.expm1(delta))
    expected = math.exp(mean+variance/2)*((1-p)+p*math.exp(delta))
    error = abs(expected-math.exp(null_variance/2))
    allowance = 128*np.finfo(float).eps*max(1., abs(expected))
    assert error <= allowance
    return dict(independent_all_allocations_minimum='4', sorted_cut_minimum=result['objective'],
        permutation_translation_scale=True, mean_rate_anchor_error=error,
        mean_rate_anchor_allowance=allowance, initializer_not_observation_model=True)


def initialize(D, ED, index, null):
    z = np.log((D.sum(axis=1)+.5)/(ED.sum(axis=1)+.5))
    centres = np.array([math.fsum(map(float, z[index == k]))/np.sum(index == k)
                        for k in range(len(null)-1)])
    record = partition(z-centres[index])
    delta, variance, p = (record[k] for k in ('delta', 'variance', 'proportion'))
    shift = np.logaddexp(math.log1p(-p), math.log(p)+delta)
    means = np.asarray(null[:-1])+.5*(null[-1]-variance)-shift
    original = np.r_[means, delta, variance]
    mirrored = np.r_[means+delta, -delta, variance]
    return (original, mirrored), dict(record, depth_centres=centres.tolist(),
        interpretation='Two-mean approximation supplies starts only; exact existing NB/state objective selects fits')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    old_task = args.task.parent/'phase8_rep2_trace_20261110'
    old = old_task/'trace_v1/replay/seed_202611098002/pooled_prevalence'
    folder = old/'guards/proposal_00000'
    assert json.loads((old_task/'summary_v1/report.json').read_text())['exact_replay_gate_passed']
    prior = json.loads((folder/'DISCOVERY_MODELS_LOCKED.json').read_text())
    assert [m['order'] for m in prior['models']] == [24, 48, 96]
    args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        exact_gate=old_task/'summary_v1/report.json', discovery=folder/'discovery.npz',
        previous_models=folder/'DISCOVERY_MODELS_LOCKED.json', previous_guard=folder/'report.json',
        fold=old/'folds/model_0.npz', objective=Path(__file__).with_name('profiled_prevalence.py'))
    for order in (24, 48, 96):
        inputs[f'prediction_{order}'] = folder/f'prediction_{order}.npz'
    run = start_run(args.out, inputs, dict(seed=202611098002, proposal=0, new_starts=2),
        dict(scope='One exposed discovery-optimizer mechanism test; no recovery, error-rate or biological validation'),
        mode='discovery_partition_start')
    try:
        checks = numerical_checks()
        write_json(checks, args.out/'numerics.json')
        with np.load(inputs['discovery']) as z:
            discovery = {key: z[key] for key in z.files}
        D, ED, rD, strata = (discovery[k] for k in ('D', 'ED', 'rD', 'strata'))
        levels, index = np.unique(strata, return_inverse=True)
        method = load('profiled_prevalence')
        models, records = [], []
        for old_model in prior['models']:
            order = old_model['order']
            assert levels.tolist() == old_model['levels']
            x = np.r_[old_model['base_mean'], old_model['delta'], old_model['variance']]
            old_value, _ = method.objective(x, D, ED, rD, index, order, shared=True)
            assert old_value == old_model['alternative_objective']
            starts, initialization = initialize(D, ED, index, old_model['null'])
            def objective(parameters):
                return method.objective(parameters, D, ED, rD, index, order, shared=True)
            fits = [minimize(objective, start, jac=True, method='L-BFGS-B',
                bounds=[(None, None)]*(len(x)-1)+[(0., None)],
                options=dict(ftol=np.finfo(float).eps, gtol=np.sqrt(np.finfo(float).eps),
                             maxiter=1000, maxls=50)) for start in starts]
            candidates = [f for f in fits if f.success and np.isfinite(f.fun)
                          and np.isfinite(f.x).all()]
            operations = 512*D.size*order
            unit = operations*np.finfo(float).eps
            assert unit < 1.
            budget = unit/(1-unit)*max(1., abs(old_value), *(abs(f.fun) for f in candidates))
            best = min(candidates, key=lambda f: f.fun) if candidates else None
            selected = bool(best is not None and best.fun+2*budget < old_value)
            model = dict(old_model)
            if selected:
                value, gradient, mixture = method.objective(best.x, D, ED, rD, index, order,
                    shared=True, details=True)
                model.update(base_mean=best.x[:-2].tolist(), delta=float(best.x[-2]),
                    variance=float(best.x[-1]), prevalence=[r['prevalence'] for r in mixture],
                    profiled_prevalence=mixture, alternative_objective=value,
                    initialization='Discovery-only exact sorted two-mean partition',
                    maximum_absolute_gradient=float(np.max(abs(gradient))))
            models.append(model)
            record = dict(order=order, selected=selected, initialization=initialization,
                old_objective=old_value, new_objective=model['alternative_objective'],
                numerical_budget=budget, attempts=[dict(success=bool(f.success),
                    objective=float(f.fun), x=f.x.tolist(), iterations=int(f.nit),
                    message=str(f.message)) for f in fits])
            records.append(record)
            write_json(dict(models=models, records=records, confirmation_counts_opened=False),
                       args.out/'DISCOVERY_MODELS_LOCKED.json')
            print(json.dumps(record), flush=True)
        # The confirmation matrix is reconstructed only AFTER every new model
        # and selection decision has been frozen from the discovery observations.
        with np.load(inputs['fold']) as z:
            def matrix(prefix):
                return sparse.csr_matrix((z[prefix+'_data'], z[prefix+'_indices'], z[prefix+'_indptr']),
                                         shape=tuple(z[prefix+'_shape']))
            positions = pd.Index(z['global_rows']).get_indexer(discovery['global_rows'])
            assert np.all(positions >= 0) and len(np.unique(positions)) == len(positions)
            columns = discovery['genes']
            source, prop = matrix('source'), matrix('prop')
            np.testing.assert_array_equal(prop[positions][:, columns].toarray(), D)
            C = (source-prop)[positions][:, columns].toarray()
            np.testing.assert_array_equal(D+C, source[positions][:, columns].toarray())
        scores, old_scores = [], []
        for old_model, model in zip(prior['models'], models):
            arguments = (D, C, ED, discovery['EC'], rD, discovery['rC'], strata)
            baseline = method.predict(*arguments, old_model)
            with np.load(inputs[f"prediction_{model['order']}"]) as z:
                for key in ('per_cell', 'discovery_membership', 'joint_membership'):
                    np.testing.assert_array_equal(baseline[key], z[key])
            previous = json.loads(inputs['previous_guard'].read_text())
            expected = next(r for r in previous['quadrature'] if r['order'] == model['order'])
            assert baseline['score'] == expected['score']
            candidate = method.predict(*arguments, model)
            arrays = {key: candidate.pop(key) for key in ('per_cell', 'discovery_membership', 'joint_membership')}
            np.savez_compressed(args.out/f"prediction_{model['order']}.npz", **arrays)
            old_scores.append(dict(order=model['order'], score=baseline['score']))
            scores.append(dict(order=model['order'], **candidate))
        diagnostic = 2*max(abs(scores[i]['score']-scores[i-1]['score']) for i in (1, 2))+sum(
            r['rounding_budget'] for r in scores)
        resolved = bool(abs(scores[-1]['score']) > diagnostic)
        passed = bool(resolved and scores[-1]['score'] > 0)
        write_json(dict(numerics=checks, fits=records, old_scores=old_scores,
            new_scores=scores, numerical_uncertainty_diagnostic=diagnostic,
            resolved=resolved, unchanged_guard_passed=passed,
            model_selection_reads_confirmation=False, old_predictions_exact=9,
            biological_specimens=0, independent_evaluation=False, original_criteria_changed=False,
            decision='Permit separately frozen recovery/control test' if passed else
                     'Reject this initializer as a sufficient correction;no seed/threshold tuning'), args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU existing profiled NB/state objective with discovery-only starts')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
