"""Independent normal-fold selection, raw-count replay, metrics and leakage audit."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from haplodose.provenance import start_run, finish_run

METHODS = ('NB_zero', 'NB_anchor', 'NB_shrink', 'pooled_zero', 'pooled_anchor', 'pooled_shrink')
CONDITIONS = ('null', 'gain', 'shared', 'mismatched', 'state')


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def same_arrays(first, second):
    a, b = np.load(first), np.load(second)
    assert set(a.files) == set(b.files)
    for key in a.files:
        assert a[key].dtype == b[key].dtype and a[key].shape == b[key].shape
        assert a[key].tobytes() == b[key].tobytes(), key
    return len(a.files)


def independent_information(data, scalar):
    """Explicit Fisher matrices, independent of the primary trace simplification."""
    probabilities = scalar.probability(data['base'], data['nuisance'], np.zeros(len(data['base'])))
    values = []
    for exposure, p in zip(data['L'], probabilities):
        mean = exposure*p
        jacobian = np.eye(len(p))-p[None, :]
        weights = (mean+data['phi']*mean**2)/(1+data['phi']*mean)**2
        fisher = jacobian.T@(weights[:, None]*jacobian)
        values.append(float(np.trace(fisher))/(len(p)-1))
    halves = [math.fsum(values[parity::2]) for parity in (0, 1)]
    return halves, math.fsum(values)


def independent_selection(data, reference, profiles, selected, information=None):
    folds = [np.arange(parity, len(data['Y']), 2) for parity in (0, 1)]
    subsets = [{k: v[rows] if k in ('Y', 'L', 'R', 'nuisance') else v for k, v in data.items()} for rows in folds]
    q0 = sum(reference.full_q(sub, sub['offset']) for sub in subsets)

    def score(t):
        def fraction(training):
            if information is None:
                return t
            mean_information = math.fsum(information)/2
            return information[training]*t/(information[training]*t+mean_information*(1-t))
        return sum(reference.full_q(subsets[held], fraction(1-held)*profiles[1-held]) for held in (0, 1))-q0

    refinements = []
    for count in (33, 65):
        grid = np.linspace(0., 1., count)
        scores = [score(t) for t in grid]
        i = int(np.argmax(scores))
        candidates = [(0., scores[0]), (1., scores[-1])]
        # A boundary grid winner can hide an interior optimum in the adjacent
        # interval; refine that interval rather than treating the mesh as truth.
        result = minimize_scalar(lambda t: -score(t), bounds=(grid[max(0, i-1)], grid[min(count-1, i+1)]),
            method='bounded', options=dict(xatol=np.finfo(float).eps**.75, maxiter=500))
        assert result.success
        candidates.append((float(result.x), score(result.x)))
        refinements.append(max(candidates, key=lambda v: v[1]))
    t, value = refinements[-1]
    selected_parameter = selected['t'] if information is None else selected['cv_parameter']
    actual = score(selected_parameter)
    eps, operations = np.finfo(float).eps, 128*data['Y'].size
    budget = operations*eps/(1-operations*eps)*max(1., abs(q0)+abs(q0+actual))
    assert abs(actual-selected['heldout_gain']) <= budget+selected['numerical_budget']
    assert abs(score(0.)-selected['endpoint_zero_gain']) <= budget
    assert abs(score(1.)-selected['endpoint_one_gain']) <= budget
    assert actual >= value-budget
    assert abs(refinements[0][1]-refinements[1][1]) <= budget
    full_fraction = t if information is None else 2*t/(1+t)
    if information is not None:
        assert abs(selected['t']-2*selected_parameter/(1+selected_parameter)) <= 128*eps
    return dict(selected_t=selected['t'], independent_t=full_fraction, independent_parameter=t, independent_gain=value,
        selected_independent_gain=actual, objective_budget=budget,
        independent_mesh_t_difference=abs(refinements[0][0]-refinements[1][0]))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--seed-key', type=int, default=904)
    p.add_argument('--information-scaled', action='store_true')
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    scalar_path = a.task.parent/'reference_offset_20261109/audit_source_v1/benchmarks/hardening/reference_offset_audit.py'
    scalar = load('scalar_baseline', scalar_path)
    reference = load('NB_reference', root/'offset_nb_numerics.py')
    old_auditor_path = a.task/'audit_source_v1/benchmarks/hardening/normal_offset_audit.py'
    old_auditor = None if a.information_scaled else load('archived_mesh_audit', old_auditor_path)
    methods = METHODS+('NB_information', 'pooled_information') if a.information_scaled else METHODS
    panel, fit = a.task/'panel_v1', a.task/'fit_v1'
    inputs = dict(script=Path(__file__), scalar=scalar_path, NB_reference=root/'offset_nb_numerics.py',
        design=a.task/'DESIGN.md', observations=panel/'observations.npz', targets=panel/'targets.npz',
        references=fit/'references.json', normal_arrays=fit/'reference_arrays.npz', fitted=fit/'fitted.npz',
        records=fit/'records.json', summary=fit/'report.json')
    if old_auditor is not None:
        inputs['archived_auditor'] = old_auditor_path
    rec = start_run(a.out, inputs, dict(seed_key=a.seed_key, information_scaled=a.information_scaled),
        dict(scope='Independent frozen conditional-RNA experiment, no biological guarantee'),
        mode='normal_offset_audit')
    try:
        leakage = []
        for variant in ('query_changed_v1', 'truth_changed_v1'):
            other = a.task/variant
            assert inputs['references'].read_bytes() == (other/'references.json').read_bytes()
            arrays = same_arrays(inputs['normal_arrays'], other/'reference_arrays.npz')
            leakage.append(dict(variant=variant, normal_arrays_bit_identical=arrays, normal_JSON_byte_identical=True))
        truth_fit_arrays = same_arrays(inputs['fitted'], a.task/'truth_changed_v1/fitted.npz')
        observed, targets = np.load(inputs['observations']), np.load(inputs['targets'])
        fitted, normal_arrays = np.load(inputs['fitted']), np.load(inputs['normal_arrays'])
        query_changed = np.load(a.task/'query_changed_v1/fitted.npz')
        assert any(fitted[k].tobytes() != query_changed[k].tobytes() for k in fitted.files)
        assert (fit/'records.json').read_bytes() != (a.task/'truth_changed_v1/records.json').read_bytes()
        references, records, summary = (json.loads(inputs[k].read_text()) for k in ('references', 'records', 'summary'))
        assert len(references) == 80 and len(records) == 80*len(methods) and len(summary['summary']) == 5*len(methods)
        selection_cache, selections, audited, refined, archived_failures = {}, [], [], [], []
        count_entries, baseline_checks = 0, 0
        eps = np.finfo(float).eps
        for rep in range(16):
            rng = np.random.default_rng(np.random.SeedSequence([20261109, a.seed_key, rep]))
            base = rng.lognormal(0, .5, 40)
            base /= base.sum()
            Ln, Lq = rng.lognormal(np.log(1000), .4, 60), rng.lognormal(np.log(2000), .4, 240)
            gamma, phi = .25*np.sin(np.arange(40)/4), np.full(40, .2)
            nn, nq = (np.log(Ln)-np.log(1000))[:, None]*gamma, (np.log(Lq)-np.log(1000))[:, None]*gamma
            r = np.zeros(40)
            r[12:22] = 1.
            gain, technical = .6*np.log(1.5)*r, .35*r

            def draw(L, nuisance, shift):
                mean = L[:, None]*scalar.probability(base, nuisance, shift)
                return rng.negative_binomial(1/phi, 1/(1+phi*mean)).astype(np.int64)

            N0, Ns = draw(Ln, nn, np.zeros(40)), draw(Ln, nn, technical)
            Q0, Qd, Qs = draw(Lq, nq, np.zeros(40)), draw(Lq, nq, gain), draw(Lq, nq, gain+technical)
            cases = dict(null=(N0, Q0, np.zeros(40), np.zeros(40), np.zeros(40)),
                gain=(N0, Qd, gain, np.zeros(40), np.zeros(40)), shared=(Ns, Qs, gain, technical, technical),
                mismatched=(Ns, Qd, gain, np.zeros(40), technical), state=(N0, Qd, np.zeros(40), gain, np.zeros(40)))
            for condition, (N, Q, dna, qother, nother) in cases.items():
                key = f'{condition}__{rep}'
                for name, value in dict(N=N, Q=Q, Ln=Ln, Lq=Lq, base=base, nn=nn, nq=nq, phi=phi).items():
                    assert np.array_equal(value, observed[key+'__'+name]), (key, name)
                count_entries += N.size+Q.size
                normal_digest = hashlib.sha256()
                for value in (N, Ln, base, nn, phi):
                    value = np.ascontiguousarray(value)
                    normal_digest.update(str((value.shape, value.dtype.str)).encode())
                    normal_digest.update(value.tobytes())
                assert references[key]['normal_input_sha256'] == normal_digest.hexdigest()
                expected_targets = dict(DNA_linked_RNA=scalar.centred(dna), RNA_zero=scalar.centred(dna+qother),
                    RNA_anchor=scalar.centred(dna+qother-nother))
                for name, value in expected_targets.items():
                    assert np.array_equal(targets[key+'__'+name], value)
                for label, rows in (('full', np.arange(60)), ('even', np.arange(0, 60, 2)), ('odd', np.arange(1, 60, 2))):
                    expected = scalar.pooled(N[rows], Ln[rows], base, nn[rows])
                    assert np.max(np.abs(expected-normal_arrays[key+'__pooled_'+label])) <= 128*N.size*eps
                    baseline_checks += 1
                data = dict(Y=N.astype(float), L=Ln, nuisance=nn, base=base, phi=phi,
                    R=np.ones((60, 1)), doses=np.zeros((1, 40)), offset=np.zeros(40))
                information = None
                if a.information_scaled:
                    information, full_information = independent_information(data, scalar)
                    information_budget = 128*N.size*eps*max(1., full_information)
                    assert abs(math.fsum(information)-full_information) <= information_budget
                    for family in ('NB', 'pooled'):
                        saved = references[key]['information'][family]
                        assert np.max(np.abs(np.array(information)-saved['fold_information'])) <= information_budget
                        assert abs(full_information-saved['full_information']) <= information_budget
                        assert abs(math.fsum(information)/2-saved['mean_fold_information']) <= information_budget
                        expected_fractions = [v*saved['cv_parameter']/(v*saved['cv_parameter']+
                            math.fsum(information)/2*(1-saved['cv_parameter'])) for v in information]
                        assert np.max(np.abs(np.array(expected_fractions)-saved['training_fractions'])) <= 128*N.size*eps
                corrected = {}
                for family, mode in [(family, mode) for family in ('NB', 'pooled')
                                     for mode in (('shrink', 'information') if a.information_scaled else ('shrink',))]:
                    normal = [normal_arrays[key+'__'+family+'_'+v] for v in ('even', 'odd')]
                    selected = references[key][family] if mode == 'shrink' else references[key]['information'][family]
                    cache_key = (references[key]['normal_input_sha256'], family, mode)
                    if cache_key not in selection_cache:
                        try:
                            if old_auditor is not None:
                                old_auditor.independent_selection(data, reference, normal, selected)
                        except AssertionError as error:
                            frame = error.__traceback__
                            while frame.tb_next is not None:
                                frame = frame.tb_next
                            values = frame.tb_frame.f_locals
                            archived_failures.append(dict(key=key, family=family,
                                selected_t=references[key][family]['t'],
                                archived_refinements=values.get('refinements'),
                                archived_budget=values.get('budget')))
                        selection_cache[cache_key] = independent_selection(data, reference, normal, selected,
                            information if mode == 'information' else None)
                    check = selection_cache[cache_key]
                    assert check['selected_t'] == selected['t']
                    selections.append(dict(key=key, family=family, mode=mode, **check))
                    query_profile = fitted[key+'__'+family+'_zero']
                    full = normal_arrays[key+'__'+family+'_full']
                    actual = scalar.centred(query_profile-selected['t']*full)
                    assert np.array_equal(actual, fitted[key+'__'+family+'_'+mode])
                    corrected[family+'_'+mode] = scalar.centred(query_profile-check['independent_t']*full)
                for method in methods:
                    estimate = fitted[key+'__'+method]
                    dna_mse = math.fsum((float(x)-float(y))**2 for x, y in zip(estimate, scalar.centred(dna)))/40
                    rna = expected_targets['RNA_zero' if method.endswith('zero') else 'RNA_anchor']
                    rna_mse = math.fsum((float(x)-float(y))**2 for x, y in zip(estimate, rna))/40
                    contrast = math.fsum(estimate[12:22])/10-math.fsum(np.r_[estimate[:12], estimate[22:]])/30
                    row = next(v for v in records if (v['condition'], v['replicate'], v['method']) == (condition, rep, method))
                    budget = 128*40*eps*max(1., dna_mse, rna_mse, abs(contrast))
                    for name, value in (('DNA_linked_RNA_MSE', dna_mse), ('own_relative_RNA_MSE', rna_mse), ('region_contrast', contrast)):
                        assert abs(row[name]-value) <= budget
                    truth_contrast = math.fsum(dna[12:22])/10-math.fsum(np.r_[dna[:12], dna[22:]])/30
                    if a.information_scaled:
                        assert abs(row['regional_squared_error']-(contrast-truth_contrast)**2) <= budget
                    audited.append(row)
                    profile = corrected.get(method, estimate)
                    refined.append(dict(condition=condition, method=method,
                        mse=float(np.mean((profile-scalar.centred(dna))**2)),
                        regional_squared_error=float((profile[12:22].mean()-np.r_[profile[:12], profile[22:]].mean()-truth_contrast)**2)))
            for method in methods:
                assert np.array_equal(fitted[f'gain__{rep}__{method}'], fitted[f'state__{rep}__{method}'])
        refined_summary = []
        for row in summary['summary']:
            values = [v for v in audited if v['condition'] == row['condition'] and v['method'] == row['method']]
            metrics = ('DNA_linked_RNA_MSE', 'own_relative_RNA_MSE', 'region_contrast')
            for name in metrics+(('regional_squared_error',) if a.information_scaled else ()):
                expected = math.fsum(v[name] for v in values)/16
                assert abs(expected-row['mean_'+name]) <= 128*eps*max(1., abs(expected))
            assert row['positive_region_contrast'] == sum(v['region_contrast'] > 0 for v in values)
            mse = math.fsum(v['mse'] for v in refined if v['condition'] == row['condition'] and v['method'] == row['method'])/16
            region_mse = math.fsum(v['regional_squared_error'] for v in refined
                if v['condition'] == row['condition'] and v['method'] == row['method'])/16
            refined_summary.append(dict(condition=row['condition'], method=row['method'], mean_DNA_linked_RNA_MSE=mse,
                mean_regional_squared_error=region_mse))
        comparisons = []
        for condition in CONDITIONS:
            pairs = [('NB_shrink', 'NB_anchor'), ('pooled_shrink', 'pooled_anchor'), ('NB_shrink', 'pooled_shrink')]
            if a.information_scaled:
                pairs += [('NB_information', 'NB_shrink'), ('pooled_information', 'pooled_shrink'),
                    ('NB_information', 'NB_anchor'), ('pooled_information', 'pooled_anchor'), ('NB_information', 'pooled_information')]
            for metric in (('mean_DNA_linked_RNA_MSE', 'mean_regional_squared_error') if a.information_scaled
                           else ('mean_DNA_linked_RNA_MSE',)):
                for first, second in pairs:
                    original = {v['method']: v[metric] for v in summary['summary'] if v['condition'] == condition}
                    alternative = {v['method']: v[metric] for v in refined_summary if v['condition'] == condition}
                    stable = (original[first] < original[second]) == (alternative[first] < alternative[second])
                    comparisons.append(dict(condition=condition, metric=metric, first=first, second=second,
                        mean_difference=original[first]-original[second], independently_refined_difference=alternative[first]-alternative[second],
                        ordering_stable=bool(stable)))
        resolved = all(v['ordering_stable'] for v in comparisons)
        (a.out/'report.json').write_text(json.dumps(dict(all_checks_passed=resolved, leakage=leakage,
            truth_changed_fit_arrays_bit_identical=truth_fit_arrays, raw_count_entries_replayed=count_entries,
            scalar_pooled_reference_checks=baseline_checks, case_method_metrics=len(audited), selection_checks=selections,
            unique_selection_inputs=len(selection_cache), summary=summary['summary'], refined_summary=refined_summary,
            archived_auditor_failures=archived_failures,
            comparisons=comparisons, biological_specimens=0, production_changed=False, scientific_criteria_changed=False,
            scope='Conditional RNA model and exposed procedure; no biological, clone or error-control claim'), indent=2)+'\n')
        if not resolved:
            raise ArithmeticError('Scientific comparison ordering is not numerically resolved')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU independent count replay, full NB grid/refinement, metrics and leakage')


if __name__ == '__main__':
    main()
