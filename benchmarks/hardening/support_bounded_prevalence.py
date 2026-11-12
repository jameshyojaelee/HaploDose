"""Probe supported-population prevalence without relaxing confirmation/support."""
import argparse
from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import minimize, minimize_scalar

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def load(name):
    spec = importlib.util.spec_from_file_location('support_bound_'+name,
        Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def interval(minimum, cells):
    if not isinstance(minimum, int) or not 0 < 2*minimum < cells:
        raise ValueError('Two supported populations require more than twice the minimum')
    exact = Fraction(minimum, cells)
    lo, hi = float(exact), float(1-exact)
    if Fraction(lo) < exact:
        lo = np.nextafter(lo, np.inf)
    if Fraction(hi) > 1-exact:
        hi = np.nextafter(hi, -np.inf)
    return lo, hi


def bounded(original, lo, hi):
    """Concavity makes clipping the exact unbounded optimum an interval optimum."""
    def fit(log0, log1):
        p, record = original(log0, log1)
        chosen = min(hi, max(lo, p))
        return chosen, dict(record, unconstrained_prevalence=p,
            population_constraint_active=bool(chosen != p),
            population_prevalence_interval=[lo, hi],
            original_score_is_unconstrained=True)
    return fit


def checks():
    method = load('profiled_prevalence')
    old = method.checks(shared=True)
    lo, hi = interval(1, 4)
    a, b = np.array([.75, .25, .75]), np.array([.25, .75, .25])
    function = bounded(method.prevalence, lo, hi)
    p, record = function(np.log(a), np.log(b))
    assert p == float(Fraction(1, 4)) and record['population_constraint_active']
    mirror, _ = function(np.log(b), np.log(a))
    assert mirror == float(Fraction(3, 4))
    flat, info = function(np.log(a), np.log(a))
    assert flat == .5 and info['identified'] is False
    # Independent full-probability scalar optimization, including both endpoints.
    objective = lambda value: -math.fsum(map(float, np.log((1-value)*a+value*b)))
    reference = minimize_scalar(objective, bounds=(lo, hi), method='bounded')
    chosen = min((lo, hi, reference.x), key=objective)
    assert chosen == p
    low, high = interval(30, 600)
    assert Fraction(low) >= Fraction(1, 20) and Fraction(high) <= Fraction(19, 20)
    # At an active bound prevalence is fixed. Compare the profiled gradient to
    # the separately verified full NB variance objective with that fixed weight.
    method.prevalence = bounded(method.prevalence, .5, .5)
    D = np.array([[0, 0], [4, 6], [0, 1], [5, 7]])
    E, r = np.tile([1., 2.], (4, 1)), np.full(D.shape, 2.)
    strata, x = np.array([0, 0, 1, 1]), np.array([-.4, -.3, .7, .04])
    value, gradient = method.objective(x, D, E, r, strata, 96, shared=True)
    ref_value, ref_gradient = 0., np.zeros_like(x)
    guard = method.load('nb_state_guard')
    for s in (0, 1):
        v, g = guard.variance_objective([x[s], x[s]+x[-2], 0., x[-1]],
            D[strata == s], E[strata == s], r[strata == s], 96, True)
        ref_value += v/2
        ref_gradient[s] = (g[0]+g[1])/2
        ref_gradient[-2:] += g[[1, 3]]/2
    assert value == ref_value
    np.testing.assert_array_equal(gradient, ref_gradient)
    return dict(existing_independent_PMF_and_gradient_checks=old,
        exact_constrained_fraction_optimum='1/4', independent_full_probability_optimum=chosen,
        exact_component_exchange=True, flat_unidentified=True,
        active_bound_objective_and_gradient_exact=True, population_bounds_round_inward=True)


def fit(method, data, old):
    D, E, r, strata = (data[k] for k in ('D', 'ED', 'rD', 'strata'))
    levels, index = np.unique(strata, return_inverse=True)
    assert levels.tolist() == old['levels']
    order, k = old['order'], len(levels)
    anchor = np.r_[old['base_mean'], old['delta'], old['variance']]
    null = np.asarray(old['null'])
    q = np.asarray(data['q'])
    pi = np.array([q[index == s].mean() for s in range(k)])
    z = np.log((D.sum(axis=1)+.5)/(E.sum(axis=1)+.5))
    spread = float(np.diff(np.quantile(z, [.2, .8]))[0])
    starts = [anchor]+[np.r_[null[:-1]-np.logaddexp(np.log1p(-pi), np.log(pi)+d), d, v]
        for d in (0., spread, -spread) for v in (0., null[-1])]
    function = lambda x: method.objective(x, D, E, r, index, order, shared=True)
    anchor_value = function(anchor)[0]
    options = dict(ftol=np.finfo(float).eps, gtol=np.sqrt(np.finfo(float).eps),
                   maxiter=1000, maxls=50)
    fits = [minimize(function, x, jac=True, method='L-BFGS-B',
        bounds=[(None, None)]*(k+1)+[(0., None)], options=options) for x in starts]
    good = [f for f in fits if f.success and np.isfinite(f.fun) and np.isfinite(f.x).all()]
    if not good:
        raise ArithmeticError('Population-constrained likelihood fit unresolved')
    best = min(good, key=lambda f: f.fun)
    unit = 512*D.size*order*np.finfo(float).eps
    assert unit < 1
    budget = unit/(1-unit)*max(1., abs(anchor_value), abs(best.fun))
    if best.fun > anchor_value+budget:
        raise ArithmeticError('Optimizer loses feasible constrained anchor likelihood')
    value, gradient, records = method.objective(best.x, D, E, r, index, order,
                                               shared=True, details=True)
    projected = gradient.copy()
    if best.x[-1] == 0:
        projected[-1] = min(0., projected[-1])
    return dict(order=order, levels=old['levels'], null=old['null'],
        base_mean=best.x[:-2].tolist(), delta=float(best.x[-2]),
        variance=float(best.x[-1]), prevalence=[rec['prevalence'] for rec in records],
        profiled_prevalence=records, alternative_objective=value,
        previous_unconstrained_objective=old['alternative_objective'],
        prevalence_shared_across_depth=True, null_unchanged=True,
        constrained_anchor_objective=anchor_value, arithmetic_budget=budget,
        maximum_projected_gradient=float(np.max(abs(projected))),
        starts=[dict(success=bool(f.success), objective=float(f.fun), x=f.x.tolist(),
            iterations=int(f.nit), message=str(f.message)) for f in fits],
        optimization='Same seven-start local fit; population prevalence constrained;global optimum not certified')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    prior = args.task.parent/'long_event_guard_probe_20261110/probe_v1'
    report = json.loads((prior/'report.json').read_text())
    assert json.loads((prior/'run_manifest.json').read_text())['status'] == 'complete'
    assert report['exact_output_parity'] == dict(arrays=4, tables=28, folds=1, exact=True)
    assert report['long_guard_passed'] is False
    short = args.task.parent/'phase8_rep2_trace_20261110/trace_v1/replay/seed_202611098002/pooled_prevalence/guards/proposal_00000'
    folders = dict(short=short, long=prior/'long')
    minimum = json.loads((folders['long']/'parameters.json').read_text())['minimum_cells']
    assert minimum == 30
    args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        prior_report=prior/'report.json', prior_manifest=prior/'run_manifest.json',
        population_parameters=folders['long']/'parameters.json',
        fold=args.task.parent/'phase8_rep2_trace_20261110/trace_v1/replay/seed_202611098002/pooled_prevalence/folds/model_0.npz',
        objective=Path(__file__).with_name('profiled_prevalence.py'))
    for name, folder in folders.items():
        for key in ('discovery.npz', 'DISCOVERY_MODELS_LOCKED.json'):
            inputs[name+'_'+key] = folder/key
        for order in (24, 48, 96):
            inputs[f'{name}_prediction_{order}'] = folder/f'prediction_{order}.npz'
    inputs['short_report'] = short/'report.json'
    run = start_run(args.out, inputs, dict(seed=202611098002, minimum=minimum,
        primary_event='event_chrom12_bins210-228', regions=['short', 'long']),
        dict(scope='Exposed constrained-population mechanism;not recovery,error control or biological validation'),
        mode='support_bounded_prevalence')
    try:
        write_json(checks(), args.out/'numerics.json')
        prepared, locked = {}, {}
        for name, folder in folders.items():
            with np.load(folder/'discovery.npz') as z:
                data = {key: z[key] for key in z.files}
            old_models = json.loads((folder/'DISCOVERY_MODELS_LOCKED.json').read_text())['models']
            assert [m['order'] for m in old_models] == [24, 48, 96]
            method = load('profiled_prevalence')
            _, index = np.unique(data['strata'], return_inverse=True)
            for model in old_models:
                actual = method.objective(np.r_[model['base_mean'], model['delta'], model['variance']],
                    data['D'], data['ED'], data['rD'], index, model['order'], shared=True)[0]
                assert actual == model['alternative_objective']
            lo, hi = interval(minimum, len(data['D']))
            method.prevalence = bounded(method.prevalence, lo, hi)
            models = []
            for old_model in old_models:
                print(json.dumps(dict(region=name, order=old_model['order'], status='discovery_fit')), flush=True)
                models.append(fit(method, data, old_model))
                locked[name] = dict(models=models, prevalence_interval=[lo, hi],
                    exact_population_floor=str(Fraction(minimum, len(data['D']))))
                write_json(dict(regions=locked, confirmation_counts_opened=False),
                    args.out/'DISCOVERY_MODELS_LOCKED.json')
            prepared[name] = data, old_models, models
        assert set(locked) == {'short', 'long'} and all(len(v['models']) == 3 for v in locked.values())
        # Only AFTER all six constrained models are locked, reconstruct C.
        with np.load(inputs['fold']) as z:
            def matrix(prefix):
                return sparse.csr_matrix((z[prefix+'_data'], z[prefix+'_indices'], z[prefix+'_indptr']),
                    shape=tuple(z[prefix+'_shape']))
            source, proposal, rows = matrix('source'), matrix('prop'), pd.Index(z['global_rows'])
        comparisons = {}
        predictor = load('profiled_prevalence')
        for name, (data, old_models, models) in prepared.items():
            positions = rows.get_indexer(data['global_rows'])
            assert rows.is_unique and np.all(positions >= 0) and len(np.unique(positions)) == len(positions)
            genes = data['genes']
            np.testing.assert_array_equal(proposal[positions][:, genes].toarray(), data['D'])
            full = source[positions][:, genes].toarray()
            C = full-data['D']
            assert np.all(C >= 0) and np.all(C == np.floor(C))
            np.testing.assert_array_equal(data['D']+C, full)
            old_scores = (json.loads((short/'report.json').read_text())['quadrature'] if name == 'short'
                          else report['records'][1]['scores'])
            scores = []
            for old_model, model, old_score in zip(old_models, models, old_scores):
                arguments = (data['D'], C, data['ED'], data['EC'], data['rD'], data['rC'], data['strata'])
                baseline = predictor.predict(*arguments, old_model)
                with np.load(inputs[f"{name}_prediction_{model['order']}"]) as z:
                    for key in ('per_cell', 'discovery_membership', 'joint_membership'):
                        np.testing.assert_array_equal(baseline[key], z[key])
                assert baseline['score'] == old_score['score']
                result = predictor.predict(*arguments, model)
                arrays = {key: result.pop(key) for key in ('per_cell', 'discovery_membership', 'joint_membership')}
                np.savez_compressed(args.out/f"{name}_prediction_{model['order']}.npz", **arrays)
                scores.append(dict(order=model['order'], **result))
            error = 2*max(abs(scores[i]['score']-scores[i-1]['score']) for i in (1, 2))+sum(r['rounding_budget'] for r in scores)
            resolved = bool(abs(scores[-1]['score']) > error)
            mass = scores[-1]['discovery_expected_members']
            comparisons[name] = dict(scores=scores, numerical_uncertainty_diagnostic=error,
                resolved=resolved, guard_passed=bool(resolved and scores[-1]['score'] > 0),
                discovery_responsibility_mass_supported=bool(min(mass) >= minimum),
                population_prior_minimum=minimum, responsibility_minimum_unchanged=minimum,
                prior_constraint_not_a_responsibility_support_substitute=True,
                old_objectives_and_prediction_arrays_and_scores_exact=True,
                confirmation_count_conservation_exact=True)
        primary = comparisons['long']
        passed = primary['guard_passed'] and primary['discovery_responsibility_mass_supported']
        write_json(dict(regions=comparisons, primary_long_event_gate_passed=bool(passed),
            criterion_changed=False, independent_evaluation=False, biological_specimens=0,
            inference_changed=False, promotion=False,
            decision='Justifies separately frozen full-workflow test only' if passed else
                'Reject population-prevalence constraint as a sufficient recovery mechanism'), args.out/'report.json')
        for value in run['input']['files'].values():
            assert file_record(value['path'])['sha256'] == value['sha256']
        finish_run(args.out, run, backend='CPU constrained shared-state NB mixture;original observation kernels')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
