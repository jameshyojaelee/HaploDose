"""Independent likelihood replay and unchanged recovery criteria for all three arms."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def recovery(task, arm, suffix):
    summary = json.loads((task/f'{arm}_{suffix}/panel/summary.json').read_text())
    assert len(summary['cases']) == 5
    outcomes = []
    for row in summary['cases']:
        seed = row['seed']
        folder = task/f'{arm}_{suffix}/panel'/f'seed_{seed}'
        truth = pd.read_csv(task.parent/f'recovery_20261107/recovery_panel_v1/seed_{seed}/data/cells.tsv', sep='\t')
        fit = pd.read_csv(folder/'fit/cells.tsv.gz', sep='\t')
        assert not any(c.startswith('truth_') for c in fit)
        assert fit.cell_id.is_unique and truth.cell_id.is_unique and set(fit.cell_id) == set(truth.cell_id)
        joined = fit.merge(truth[['cell_id', 'truth_clone']], on='cell_id', validate='one_to_one')
        s1 = joined[joined.specimen.eq('S1')]
        normal = s1.truth_clone.eq('normal')
        fpr = sum(s1.loc[normal, 'malignant_call'])/sum(normal)
        tpr = sum(s1.loc[~normal, 'malignant_call'])/sum(~normal)
        s2 = joined[(joined.specimen == 'S2') & joined.malignant_call]
        sub = s1[s1.truth_clone.isin(['clone1', 'clone2']) & s1.malignant_call]
        groups, populations = sorted(set(sub.clone)), sorted(set(sub.truth_clone))
        table = [[int(sum((sub.truth_clone == pop) & (sub.clone == group))) for group in groups] for pop in populations]
        n = min(len(groups), len(populations))
        diagonals = [sum(table[i][i] for i in range(n)), sum(table[len(populations)-1-i][i] for i in range(n))]
        accuracy = max(diagonals)/len(sub) if len(sub) else None
        differences = pd.read_csv(folder/'fit/differences.tsv', sep='\t')
        differences = differences[(differences.specimen == 'S1') & ~differences.clone_a.str.endswith('_c0')
                                  & ~differences.clone_b.str.endswith('_c0')]
        chromosomes = sorted(set(differences.chromosome))
        expected = dict(normal_fpr=(fpr, fpr < .1), tumour_tpr=(tpr, tpr > .9),
            s2_one_group=(len(set(s2.clone)), len(set(s2.clone)) == 1), s1_two_groups=(len(groups), len(groups) == 2),
            assignment_accuracy=(accuracy, accuracy is not None and accuracy > .65),
            difference_chromosomes=(chromosomes, chromosomes == ['chr12']))
        for name, (value, passed) in expected.items():
            assert row['criteria'][name]['value'] == value
            assert bool(row['criteria'][name]['passed']) == bool(passed)
        assert row['original_assertion_passed'] == all(v[1] for v in expected.values())
        assert json.loads((folder/'fit/qc.json').read_text())['params']['split_mode'] == 'thin'
        outcomes.append(dict(arm=arm, seed=seed, passed=row['original_assertion_passed']))
    return summary, outcomes


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--panel-suffix', default='stable_v1')
    p.add_argument('--numerical-gate', default='numerics_stable_v1')
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    reference = load('reference', root/'offset_nb_numerics.py')
    helper = load('conditional_objective', root/'offset_nb_objective.py')
    inputs = dict(script=Path(__file__), reference=root/'offset_nb_numerics.py', helper=root/'offset_nb_objective.py',
        design=a.task/'DESIGN.md', numerical_gate=a.task/a.numerical_gate/'report.json')
    for arm in ('legacy', 'conditional', 'NB'):
        inputs[f'{arm}_summary'] = a.task/f'{arm}_{a.panel_suffix}/panel/summary.json'
        inputs[f'{arm}_trace'] = a.task/f'{arm}_{a.panel_suffix}/trace.json'
        inputs[f'{arm}_contexts'] = a.task/f'{arm}_{a.panel_suffix}/contexts.json'
    rec = start_run(a.out, inputs, {}, dict(scope='Independent conditional-Q replay; all original criteria'),
        mode='offset_objective_audit')
    try:
        assert json.loads(inputs['numerical_gate'].read_text())['all_checks_passed']
        all_outcomes, audited, traces = [], [], []
        eps = np.finfo(float).eps
        for arm in ('legacy', 'conditional', 'NB'):
            summary, outcomes = recovery(a.task, arm, a.panel_suffix)
            all_outcomes.extend(outcomes)
            if arm != 'NB':
                old_name = 'recovery_population_v1' if arm == 'legacy' else 'recovery_after_v1'
                old = json.loads((a.task.parent/'offset_expectation_20261109'/old_name/'summary.json').read_text())
                assert [(v['seed'], v['criteria'], v['original_assertion_passed']) for v in summary['cases']] == [
                    (v['seed'], v['criteria'], v['original_assertion_passed']) for v in old['cases']]
            trace = json.loads(inputs[f'{arm}_trace'].read_text())
            contexts = json.loads(inputs[f'{arm}_contexts'].read_text())
            for row in trace:
                assert row['decreased'] == (row['gain'] < -row['budget'])
            for item in contexts:
                with np.load(item['path']) as saved:
                    data = {k: saved[k] for k in ('Y', 'L', 'R', 'doses', 'base', 'nuisance', 'phi', 'offset')}
                    offset, doses = saved['candidate_offset'], saved['candidate_doses']
                    q0, q1 = reference.full_q(data, data['offset']), reference.full_q(data, offset, doses)
                    gain, _, helper_budget = helper.gain_gradient(data, offset, doses)
                    ops = 128*data['Y'].size*data['R'].shape[1]
                    reference_budget = ops*eps/(1-ops*eps)*max(1., abs(q0)+abs(q1))
                    assert abs((q1-q0)-gain) <= reference_budget+helper_budget
                    assert item['row']['gain'] == gain
                    assert item['row'] == min((v for v in trace if v['model'] == item['row']['model'] and
                                              v['kind'] == item['row']['kind']), key=lambda v: v['gain'])
                    idx = np.unique([0, len(data['Y'])-1])
                    subset = {k: v[idx] if k in ('Y', 'L', 'R', 'nuisance') else v for k, v in data.items()}
                    sample_gain, _, sample_budget = helper.gain_gradient(subset, offset, doses)
                    decimal = reference.decimal_gain(subset, offset, doses)
                    assert abs(sample_gain-decimal) <= sample_budget
                    audited.append(dict(arm=arm, model=item['row']['model'], kind=item['row']['kind'],
                        full_NB_gain=q1-q0, helper_gain=gain, reference_budget=reference_budget,
                        subset_Decimal_gain=decimal, subset_error=abs(sample_gain-decimal),
                        subset_budget=sample_budget,
                        independently_negative=bool((q1-q0) < -(reference_budget+helper_budget))))
            traces.append(dict(arm=arm, updates=len(trace), contexts=len(contexts),
                profile_decreases=sum(v['decreased'] for v in trace if v['kind'] == 'profile'),
                offset_decreases=sum(v['decreased'] for v in trace if v['kind'] == 'offset'),
                optimizer_failures=sum(not v['optimizer']['success'] for v in trace if 'optimizer' in v)))
        assert len(all_outcomes) == 15
        (a.out/'report.json').write_text(json.dumps(dict(recovery=all_outcomes, independent_criteria=90,
            contexts=audited, traces=traces, all_checks_passed=True, Decimal_precision=80,
            biological_specimens=0, production_changed=False, scientific_criteria_changed=False), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU independent SciPy full NB, bounded Decimal80 and recovery recount')


if __name__ == '__main__':
    main()
