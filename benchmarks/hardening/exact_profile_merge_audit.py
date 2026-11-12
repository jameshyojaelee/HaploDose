"""Independently recount fixed mechanism endpoints and inspect exact redundancy."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def recount(fit, truth, condition):
    c = pd.read_csv(fit/'cells.tsv.gz', sep='\t')
    c = c[c.specimen.eq('S1')].merge(truth[['cell_id', 'truth_clone']], on='cell_id', validate='one_to_one')
    assert len(c) == 600 and c.cell_id.is_unique and not c.truth_clone.isna().any()
    table = pd.read_csv(fit/'clones.tsv', sep='\t')
    table = table[~table.normal_like & table.specimen.eq('S1')]
    groups = sorted(table.clone)
    pops = sorted(set(c.truth_clone))
    matrix = [[int(sum(c.truth_clone.eq(pop) & c.clone.eq(g) & c.malignant_call)) for g in groups] for pop in pops]
    n = min(len(pops), len(groups))
    correct = max((sum(matrix[i][j] for i, j in zip(rows, cols))
                   for rows in itertools.permutations(range(len(pops)), n)
                   for cols in itertools.combinations(range(len(groups)), n)), default=0)
    d = pd.read_csv(fit/'differences.tsv', sep='\t')
    chrom = sorted(set(d.loc[d.specimen.eq('S1') & d.clone_a.isin(groups) & d.clone_b.isin(groups), 'chromosome']))
    return dict(reported_groups=len(groups), null_error=len(groups)>1,
                exact_chr12_recovery=len(groups)==2 and chrom==['chr12'],
                assignment_correct=correct, assignment_n=600, tumour_calls=int(c.malignant_call.sum()),
                difference_chromosomes=chrom, assignment_table=matrix, groups=groups, populations=pops)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--replay-directory', default='mechanisms_v1')
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    source = a.task/a.replay_directory
    rows = json.loads((source/'report.json').read_text())
    prior = json.loads((a.task.parent/'gene_state_trace_20261109/trace_v1/report.json').read_text())
    rec = start_run(a.out, dict(script=Path(__file__), report=source/'report.json', design=a.task/'DESIGN.md'), {},
                    dict(scope='Exposed six-case mechanisms, not new biological error calibration'), mode='exact_merge_audit')
    try:
        result = []
        for row, old in zip(rows, prior):
            assert (row['condition'], row['rep']) == (old['condition'], old['rep'])
            truth = pd.read_csv(Path(old['prior'])/'input/truth_cells.tsv', sep='\t')
            for key, fit in [('before', Path(old['prior'])/'candidate/fit'), ('after', Path(row['out'])/'fit')]:
                for name, value in recount(fit, truth, row['condition']).items():
                    assert row[key][name] == value, (row['condition'], row['rep'], key, name)
            trace = json.loads((Path(row['out'])/'trace.json').read_text())
            final = next(e for e in trace['stages'] if e['stage']=='final')
            with np.load(final['arrays']) as z:
                dose, R = z['doses'], z['R']
                assert all(not np.array_equal(dose[i], dose[j])
                           for i, j in itertools.combinations(range(1, len(dose)), 2))
                assert R.shape[1] == len(dose)
            merges = [e for e in final['log'] if e['action']=='merge_identical_profiles']
            algebra = []
            for call in trace['merges']:
                with np.load(call['arrays']) as z:
                    old_dose = z['doses_before']
                    if call['statistical_merge']:
                        continue
                    groups = [[0]]
                    for k in range(1, len(old_dose)):
                        same = next((g for g in groups[1:] if np.isfinite(old_dose[k]).all()
                                     and np.array_equal(old_dose[k], old_dose[g[0]])), None)
                        if same is None:
                            groups.append([k])
                        else:
                            same.append(k)
                    expected_pi = np.array([z['pi_before'][g].sum() for g in groups])
                    expected_R = np.column_stack([z['R_before'][:, g].sum(axis=1) for g in groups])
                    budget = 64 * z['R_before'].size * np.finfo(z['R_before'].dtype).eps
                    np.testing.assert_array_equal(z['doses_after'], old_dose[[g[0] for g in groups]])
                    np.testing.assert_allclose(z['pi_after'], expected_pi/expected_pi.sum(), rtol=0, atol=budget)
                    np.testing.assert_allclose(z['R_after'], expected_R, rtol=0, atol=budget)
                    # Independent normalized rate mixture; row-wise shared
                    # nuisance terms cannot distinguish exact dose vectors.
                    lograte = old_dose-logsumexp(old_dose, axis=1, keepdims=True)
                    lhs = logsumexp(lograte+np.log(z['pi_before'])[:, None], axis=0)
                    rhs = logsumexp(lograte[[g[0] for g in groups]]+np.log(expected_pi)[:, None], axis=0)
                    np.testing.assert_allclose(lhs, rhs, rtol=0, atol=budget)
                    algebra.append(dict(before=len(old_dose), after=len(groups), invariant_checks_passed=True))
            result.append(dict(condition=row['condition'], rep=row['rep'], before=row['before'], after=row['after'],
                               exact_merges=merges, algebra=algebra, no_final_exact_tumour_duplicates=True))
        # This is the originally identified redundant four-component specimen.
        target = next(e for e in result if e['condition']=='state' and e['rep']==1)
        assert target['before']['reported_groups']==4 and target['after']['reported_groups']==1
        assert any(e['before']>e['after'] for e in target['algebra'])
        write_json(dict(all_checks_passed=True, cases=result, biological_specimens=0,
                        scientific_criteria_changed=False,
                        interpretation='Computational repair; state calibration and biological validity unresolved'), a.out/'report.json')
        finish_run(a.out, rec, backend='CPU exhaustive assignment and public-output recount; exact latent equality')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
