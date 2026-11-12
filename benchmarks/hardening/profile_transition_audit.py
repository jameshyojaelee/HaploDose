"""Audit applied profile transitions, exact baseline replay and original recovery."""
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
    recovery = load('recovery', root/'offset_objective_audit.py')
    fit = a.task/'legacy_v1'
    previous = a.task.parent/'offset_objective_20261109/legacy_stable_v1'
    inputs = dict(script=Path(__file__), reference=root/'offset_nb_numerics.py',
        helper=root/'offset_nb_objective.py', recovery=root/'offset_objective_audit.py',
        design=a.task/'DESIGN.md', trace=fit/'trace.json', contexts=fit/'contexts.json',
        summary=fit/'panel/summary.json', old_trace=previous/'trace.json',
        old_summary=previous/'panel/summary.json')
    rec = start_run(a.out, inputs, {}, dict(scope='Executed conditional profile transition; biological n=0'),
        mode='profile_transition_audit')
    try:
        trace = json.loads(inputs['trace'].read_text())
        old = json.loads(inputs['old_trace'].read_text())
        raw = [row for row in trace if row['kind'] != 'profile_applied']
        assert len(raw) == len(old)
        for index, (actual, expected) in enumerate(zip(raw, old)):
            assert dict(actual, index=index) == expected
        summary, outcomes = recovery.recovery(a.task, 'legacy', 'v1')
        old_summary = json.loads(inputs['old_summary'].read_text())
        assert [(v['seed'], v['criteria'], v['original_assertion_passed']) for v in summary['cases']] == [
            (v['seed'], v['criteria'], v['original_assertion_passed']) for v in old_summary['cases']]
        raw_profiles = [row for row in trace if row['kind'] == 'profile']
        applied = [row for row in trace if row['kind'] == 'profile_applied']
        assert len(raw_profiles) == len(applied)
        transitions = []
        for original, actual in zip(raw_profiles, applied):
            assert original['model'] == actual['model']
            assert original['components'] == actual['old_components']
            assert actual['components'] == len(actual['keep'])
            assert actual['keep'][0] == 0 and len(set(actual['keep'])) == len(actual['keep'])
            assert min(actual['keep']) >= 0 and max(actual['keep']) < actual['old_components']
            assert actual['actual_normal_dose_zero'] and actual['retained_mass_not_renormalized']
            transitions.append(dict(model=actual['model'], raw_gain=original['gain'], applied_gain=actual['gain'],
                raw_decreased=original['decreased'], applied_decreased=actual['decreased'],
                dropped_components=actual['old_components']-actual['components']))
        contexts = json.loads(inputs['contexts'].read_text())
        audited = []
        eps = np.finfo(float).eps
        for item in contexts:
            row = item['row']
            saved = np.load(item['path'])
            data = {k: saved[k] for k in ('Y', 'L', 'R', 'doses', 'base', 'nuisance', 'phi', 'offset')}
            candidate, doses = saved['candidate_offset'], saved['candidate_doses']
            q0, q1 = reference.full_q(data, data['offset']), reference.full_q(data, candidate, doses)
            gain, _, helper_budget = helper.gain_gradient(data, candidate, doses)
            operations = 128*data['Y'].size*data['R'].shape[1]
            budget = operations*eps/(1-operations*eps)*max(1., abs(q0)+abs(q1))
            assert abs((q1-q0)-gain) <= budget+helper_budget
            assert row['gain'] == gain
            assert row == min((v for v in trace if v['model'] == row['model'] and
                v['kind'] == row['kind']), key=lambda v: v['gain'])
            rows = np.unique([0, len(data['Y'])-1])
            subset = {k: v[rows] if k in ('Y', 'L', 'R', 'nuisance') else v for k, v in data.items()}
            decimal = reference.decimal_gain(subset, candidate, doses)
            sample_gain, _, sample_budget = helper.gain_gradient(subset, candidate, doses)
            assert abs(sample_gain-decimal) <= sample_budget
            if row['kind'] == 'profile_applied':
                assert np.array_equal(data['doses'][0], np.zeros(data['Y'].shape[1]))
                assert np.array_equal(doses[0], np.zeros(data['Y'].shape[1]))
                assert np.array_equal(saved['keep'], row['keep'])
                assert np.array_equal(data['R'], saved['pre_keep_R'][:, saved['keep']])
                assert np.array_equal(data['doses'], saved['pre_keep_doses'][saved['keep']])
                # R is returned at native precision before conversion for scoring.
                unit = np.finfo(np.dtype(row['responsibility_dtype'])).eps
                operations_R = row['old_components']*4
                responsibility_budget = operations_R*unit/(1-operations_R*unit)
                assert np.isfinite(data['R']).all() and np.all(data['R'] >= 0)
                assert np.all(data['R'].sum(axis=1) <= 1+responsibility_budget)
            audited.append(dict(model=row['model'], kind=row['kind'], primary_gain=gain,
                independent_gain=q1-q0, accumulation_budget=budget, helper_budget=helper_budget,
                Decimal80_gain=decimal, sample_gain=sample_gain, sample_budget=sample_budget))
        reductions = []
        for kind in ('profile', 'profile_applied', 'offset'):
            values = [row for row in trace if row['kind'] == kind]
            reductions.append(dict(kind=kind, updates=len(values), decreased=sum(v['decreased'] for v in values)))
        (a.out/'report.json').write_text(json.dumps(dict(all_checks_passed=True,
            original_criteria_checked=30, recovery=outcomes, raw_trace_rows_exactly_replayed=len(raw),
            summary=reductions, transitions=transitions, contexts=audited, Decimal_precision=80,
            biological_specimens=0, production_changed=False, scientific_criteria_changed=False,
            scope='Conditional NB transition diagnostic; a decrease does not by itself invalidate pooled PLN segmentation'),
            indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU independent full NB, Decimal80, trace replay and30 recovery criteria')


if __name__ == '__main__':
    main()
