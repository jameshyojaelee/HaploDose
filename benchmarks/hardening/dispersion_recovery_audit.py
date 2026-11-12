"""Independently recount all original recovery criteria and local moment updates."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.provenance import start_run, finish_run


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    rec = start_run(a.out, dict(script=Path(__file__), design=a.task/'RECOVERY_DESIGN.md',
        MAP=a.task/'recovery_MAP_v1/summary.json', soft=a.task/'recovery_soft_v1/summary.json'), {},
        dict(scope='Original criteria unchanged; independent recount, no biological claim'), mode='dispersion_recovery_audit')
    try:
        checked, update_counts = [], []
        for arm in ('MAP', 'soft'):
            summary = json.loads((a.task/f'recovery_{arm}_v1/summary.json').read_text())
            assert len(summary['cases']) == 5
            for row in summary['cases']:
                seed = row['seed']
                folder = a.task/f'recovery_{arm}_v1'/f'seed_{seed}'
                assert json.loads((folder/'fit/qc.json').read_text())['params']['split_mode'] == 'thin'
                truth = pd.read_csv(a.task.parent/f'recovery_20261107/recovery_panel_v1/seed_{seed}/data/cells.tsv', sep='\t')
                fit = pd.read_csv(folder/'fit/cells.tsv.gz', sep='\t')
                assert not any(c.startswith('truth_') for c in fit)
                assert fit.cell_id.is_unique and truth.cell_id.is_unique
                fit = fit.merge(truth[['cell_id', 'truth_clone']], on='cell_id', validate='one_to_one')
                s1 = fit[fit.specimen.eq('S1')]
                normal = s1.truth_clone.eq('normal')
                fpr = sum(bool(v) for v in s1.loc[normal, 'malignant_call'])/sum(normal)
                tpr = sum(bool(v) for v in s1.loc[~normal, 'malignant_call'])/sum(~normal)
                s2 = fit[(fit.specimen == 'S2') & fit.malignant_call]
                sub = s1[s1.truth_clone.isin(['clone1', 'clone2']) & s1.malignant_call]
                groups = sorted(set(sub.clone))
                populations = sorted(set(sub.truth_clone))
                matrix = [[int(sum((sub.truth_clone == population) & (sub.clone == group))) for group in groups]
                          for population in populations]
                forward = sum(matrix[i][i] for i in range(min(len(populations), len(groups))))
                reverse = sum(matrix[len(populations)-1-i][i] for i in range(min(len(populations), len(groups))))
                accuracy = max(forward, reverse)/len(sub) if len(sub) else None
                differences = pd.read_csv(folder/'fit/differences.tsv', sep='\t')
                differences = differences[(differences.specimen == 'S1') & ~differences.clone_a.str.endswith('_c0')
                                          & ~differences.clone_b.str.endswith('_c0')]
                chromosomes = sorted(set(differences.chromosome))
                expected = dict(normal_fpr=(fpr, fpr < .1), tumour_tpr=(tpr, tpr > .9),
                    s2_one_group=(len(set(s2.clone)), len(set(s2.clone)) == 1),
                    s1_two_groups=(len(groups), len(groups) == 2), assignment_accuracy=(accuracy, accuracy is not None and accuracy > .65),
                    difference_chromosomes=(chromosomes, chromosomes == ['chr12']))
                for name, (value, passed) in expected.items():
                    assert row['criteria'][name]['value'] == value
                    assert bool(row['criteria'][name]['passed']) == bool(passed)
                assert bool(row['original_assertion_passed']) == all(v[1] for v in expected.values())
                trace = json.loads((folder/'dispersion_trace.json').read_text())
                assert (len(trace) > 0) == (arm == 'soft')
                for update in trace:
                    num, den = np.array(update['numerator']), np.array(update['denominator'])
                    moment = np.clip(np.divide(num, den, out=np.full(len(num), .05), where=den > 0), .001, 10.)
                    phi = .5*moment+.5*np.median(moment)
                    assert np.array_equal(phi, update['soft_phi'])
                    assert np.median(phi) == update['soft_phi_median']
                checked.append(dict(arm=arm, seed=seed, passed=bool(row['original_assertion_passed']), criteria=6))
                update_counts.append(dict(arm=arm, seed=seed, updates=len(trace)))
        (a.out/'report.json').write_text(json.dumps(dict(cases=checked, update_counts=update_counts,
            independent_criteria=60, all_checks_passed=True, biological_specimens=0,
            production_changed=False, scientific_criteria_changed=False), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU independent count/identity criteria and dispersion reconstruction')


if __name__ == '__main__':
    main()
