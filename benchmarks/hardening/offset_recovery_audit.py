"""Recount the unchanged six recovery criteria from cell identities and contrasts."""
import argparse
import json
from pathlib import Path

import pandas as pd
from haplodose.provenance import start_run, finish_run


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--include-production', action='store_true')
    p.add_argument('--arms', nargs='+', help='Explicit replay directories; original criteria remain fixed')
    a = p.parse_args()
    arms = a.arms if a.arms is not None else ['recovery_before_v1', 'recovery_after_v1', 'population_only_v1/panel', 'nuisance_only_v1/panel']
    if a.include_production:
        arms.append('recovery_population_v1')
    a.out.mkdir(exist_ok=False)
    rec = start_run(a.out, dict(script=Path(__file__),
        **{arm: a.task/arm/'summary.json' for arm in arms}), {},
        dict(scope='Development-only original recovery criteria; no biological calibration'), mode='offset_recovery_audit')
    try:
        checked = []
        for arm in arms:
            summary = json.loads((a.task/arm/'summary.json').read_text())
            assert len(summary['cases']) == 5
            for row in summary['cases']:
                seed = row['seed']
                folder = a.task/arm/f'seed_{seed}'
                qc = json.loads((folder/'fit/qc.json').read_text())
                assert qc['params']['split_mode'] == 'thin'
                truth = pd.read_csv(a.task.parent/f'recovery_20261107/recovery_panel_v1/seed_{seed}/data/cells.tsv', sep='\t')
                fit = pd.read_csv(folder/'fit/cells.tsv.gz', sep='\t')
                assert not any(c.startswith('truth_') for c in fit)
                assert fit.cell_id.is_unique and truth.cell_id.is_unique
                assert set(fit.cell_id) == set(truth.cell_id)
                fit = fit.merge(truth[['cell_id', 'truth_clone']], on='cell_id', validate='one_to_one')
                s1 = fit[fit.specimen.eq('S1')]
                normal = s1.truth_clone.eq('normal')
                fpr = sum(s1.loc[normal, 'malignant_call'])/sum(normal)
                tpr = sum(s1.loc[~normal, 'malignant_call'])/sum(~normal)
                s2 = fit[(fit.specimen == 'S2') & fit.malignant_call]
                sub = s1[s1.truth_clone.isin(['clone1', 'clone2']) & s1.malignant_call]
                groups, populations = sorted(set(sub.clone)), sorted(set(sub.truth_clone))
                matrix = [[int(sum((sub.truth_clone == pop) & (sub.clone == group))) for group in groups] for pop in populations]
                n = min(len(groups), len(populations))
                diagonals = [sum(matrix[i][i] for i in range(n)),
                             sum(matrix[len(populations)-1-i][i] for i in range(n))]
                accuracy = max(diagonals)/len(sub) if len(sub) else None
                differences = pd.read_csv(folder/'fit/differences.tsv', sep='\t')
                differences = differences[(differences.specimen == 'S1') & ~differences.clone_a.str.endswith('_c0')
                                          & ~differences.clone_b.str.endswith('_c0')]
                chromosomes = sorted(set(differences.chromosome))
                expected = dict(normal_fpr=(fpr, fpr < .1), tumour_tpr=(tpr, tpr > .9),
                    s2_one_group=(len(set(s2.clone)), len(set(s2.clone)) == 1),
                    s1_two_groups=(len(groups), len(groups) == 2),
                    assignment_accuracy=(accuracy, accuracy is not None and accuracy > .65),
                    difference_chromosomes=(chromosomes, chromosomes == ['chr12']))
                for name, (value, passed) in expected.items():
                    assert row['criteria'][name]['value'] == value
                    assert bool(row['criteria'][name]['passed']) == bool(passed)
                passed = all(v[1] for v in expected.values())
                assert bool(row['original_assertion_passed']) == passed
                checked.append(dict(arm=arm, seed=seed, passed=passed, criteria=6))
        (a.out/'report.json').write_text(json.dumps(dict(cases=checked,
            independent_criteria=len(checked)*6, all_checks_passed=True,
            biological_specimens=0, scientific_criteria_changed=False), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU independent cell-identity and difference recount')


if __name__ == '__main__':
    main()
