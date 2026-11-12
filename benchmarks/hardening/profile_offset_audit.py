"""Independent profile mean, conditional-Q and unchanged recovery audit."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    reference = load('independent_NB', root/'offset_nb_numerics.py')
    helper = load('primary_NB', root/'offset_nb_objective.py')
    forward = load('independent_means', root/'profile_offset_numerics.py')
    recovery = load('recovery', root/'offset_objective_audit.py')
    inputs = dict(script=Path(__file__), reference=root/'offset_nb_numerics.py',
        helper=root/'offset_nb_objective.py', forward=root/'profile_offset_numerics.py',
        recovery=root/'offset_objective_audit.py', design=a.task/'DESIGN.md',
        gate=a.task/'numerics_v2/report.json')
    for arm in ('legacy', 'NB'):
        folder = a.task/f'profile_{arm}_v1'
        for name in ('contexts', 'trace', 'profile_trace', 'expected_contexts', 'report'):
            inputs[arm+'_'+name] = folder/(name+'.json')
        inputs[arm+'_summary'] = folder/'panel/summary.json'
        inputs[arm+'_comparator'] = a.task.parent/f'offset_objective_20261109/{arm}_stable_v1/panel/summary.json'
    rec = start_run(a.out, inputs, {}, dict(scope='Conditional forward/profile mechanism and original recovery, biological n=0'),
        mode='profile_offset_audit')
    try:
        assert json.loads(inputs['gate'].read_text())['all_checks_passed']
        eps = np.finfo(float).eps
        outcomes, comparators, means, likelihoods, reductions, profile_calls = [], [], [], [], [], []
        for arm in ('legacy', 'NB'):
            summary, independent_outcomes = recovery.recovery(a.task, 'profile_'+arm, 'v1')
            outcomes.extend(independent_outcomes)
            old = json.loads(inputs[arm+'_comparator'].read_text())
            assert [v['seed'] for v in old['cases']] == [v['seed'] for v in summary['cases']]
            comparators.extend(dict(arm=arm+'_only', seed=v['seed'], passed=v['original_assertion_passed'],
                criteria=v['criteria']) for v in old['cases'])
            trace = json.loads(inputs[arm+'_trace'].read_text())
            calls = json.loads(inputs[arm+'_profile_trace'].read_text())
            assert calls and all(v['original_base_restored'] and v['original_gamma_restored'] and
                v['offset_unchanged'] for v in calls)
            profile_calls.append(dict(arm=arm, calls=len(calls), restored=len(calls)))
            expected_contexts = json.loads(inputs[arm+'_expected_contexts'].read_text())
            assert len({(v['model'], v['nonconstant_offset']) for v in expected_contexts}) == len(expected_contexts)
            for item in expected_contexts:
                saved = np.load(item['path'])
                values = forward.independent_expected(saved['base'], saved['offset'], saved['original_gamma'],
                    saved['Lk'], saved['XLk'])
                operations = 128*values.size
                budget = operations*eps/(1-operations*eps)*max(1., float(np.max(saved['Lk'])))
                error = float(np.max(np.abs(values-saved['expected'])))
                assert error <= budget
                assert np.max(np.abs(saved['expected'].sum(axis=1)-saved['Lk'])) <= budget
                probability = forward.independent_expected(saved['base'], saved['offset'],
                    np.zeros(len(saved['base'])), np.array([1.]), np.zeros((1, 2)))[0]
                assert np.max(np.abs(probability-saved['profile_base'])) <= operations*eps
                gamma = saved['original_gamma'].copy()
                if item['nonconstant_offset']:
                    gamma -= float(probability@gamma)
                assert np.max(np.abs(gamma-saved['gamma'])) <= operations*eps*max(1., float(np.max(np.abs(gamma))))
                means.append(dict(arm=arm, model=item['model'], nonconstant_offset=item['nonconstant_offset'],
                    max_mean_error=error, budget=budget))
            applied = [v for v in trace if v['kind'] == 'profile_applied']
            assert len(applied) == sum(v['kind'] == 'profile' for v in trace)
            for row in applied:
                assert row['actual_normal_dose_zero'] and row['retained_mass_not_renormalized']
                assert row['keep'][0] == 0 and len(set(row['keep'])) == len(row['keep']) == row['components']
            for item in json.loads(inputs[arm+'_contexts'].read_text()):
                row = item['row']
                saved = np.load(item['path'])
                data = {k: saved[k] for k in ('Y', 'L', 'R', 'doses', 'base', 'nuisance', 'phi', 'offset')}
                offset, doses = saved['candidate_offset'], saved['candidate_doses']
                q0, q1 = reference.full_q(data, data['offset']), reference.full_q(data, offset, doses)
                gain, _, primary_budget = helper.gain_gradient(data, offset, doses)
                operations = 128*data['Y'].size*data['R'].shape[1]
                budget = operations*eps/(1-operations*eps)*max(1., abs(q0)+abs(q1))
                assert abs(q1-q0-gain) <= budget+primary_budget
                assert gain == row['gain']
                assert row == min((v for v in trace if v['model'] == row['model'] and v['kind'] == row['kind']),
                    key=lambda v: v['gain'])
                cells = np.unique([0, len(data['Y'])-1])
                sub = {k: v[cells] if k in ('Y', 'L', 'R', 'nuisance') else v for k, v in data.items()}
                decimal = reference.decimal_gain(sub, offset, doses)
                sample, _, sample_budget = helper.gain_gradient(sub, offset, doses)
                assert abs(sample-decimal) <= sample_budget
                if row['kind'] == 'profile_applied':
                    assert np.array_equal(data['R'], saved['pre_keep_R'][:, saved['keep']])
                    assert np.array_equal(data['doses'], saved['pre_keep_doses'][saved['keep']])
                    assert np.array_equal(doses[0], np.zeros(data['Y'].shape[1]))
                likelihoods.append(dict(arm=arm, model=row['model'], kind=row['kind'], gain=gain,
                    independent_gain=q1-q0, budget=budget, primary_budget=primary_budget,
                    subset_Decimal80_gain=decimal, subset_gain=sample, subset_budget=sample_budget))
            report = json.loads(inputs[arm+'_report'].read_text())
            for row in report['summary']:
                values = [v for v in trace if v['kind'] == row['kind']]
                assert row['updates'] == len(values) and row['decreased'] == sum(v['decreased'] for v in values)
                assert all(v['decreased'] == (v['gain'] < -v['budget']) for v in values)
                reductions.append(dict(arm=arm, **row))
        assert len(outcomes) == 10
        (a.out/'report.json').write_text(json.dumps(dict(all_checks_passed=True, original_criteria_checked=60,
            recovery=outcomes, frozen_comparators=comparators, means=means, likelihood_contexts=likelihoods,
            summary=reductions, profile_calls=profile_calls, biological_specimens=0,
            production_changed=False, scientific_criteria_changed=False,
            scope='Original development inputs and conditional model checks; no new biological/null/power validation'),
            indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU independent Decimal80 profile means/NB ratios and full SciPy NB Q/recovery')


if __name__ == '__main__':
    main()
