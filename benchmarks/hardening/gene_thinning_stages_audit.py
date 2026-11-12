"""Independent recount of unchanged recovery and recorded membership transitions."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.io import write_json
from haplodose.provenance import file_record, start_run, finish_run


SEEDS = (20261108, 20261109, 20261110)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=a.task/'DESIGN.md')
    for seed in SEEDS:
        folder = a.task/f'seed_{seed}'
        for path in sorted(folder.rglob('*')):
            if path.is_file() and path.suffix in ('.json', '.npz', '.tsv', '.gz'):
                inputs[f'{seed}:{path.relative_to(folder)}'] = path
        src = a.task.parent/f'recovery_20261107/recovery_panel_v1/seed_{seed}'
        inputs[f'{seed}:truth'] = src/'data/cells.tsv'
        inputs[f'{seed}:metadata'] = src/'model_cells.tsv'
        previous = a.task.parent/f'gene_thinning_20261109/gene_v1/panel/seed_{seed}'
        for name in ('cells.tsv.gz', 'differences.tsv'):
            inputs[f'{seed}:previous:{name}'] = previous/'fit'/name
        inputs[f'{seed}:previous:summary'] = previous/'summary.json'
    rec = start_run(a.out, inputs, {}, dict(scope='Postfit development mechanism only, biological n=0'),
        mode='gene_thinning_stages_audit')
    summaries, candidates, memberships, events_out = [], [], [], []
    try:
        for seed in SEEDS:
            folder = a.task/f'seed_{seed}'
            trace = json.loads((folder/'trace.json').read_text())
            report = json.loads((folder/'report.json').read_text())
            assert not report['errors'] and report['fit_tables_replay_exactly']
            for record in report['artifacts'].values():
                assert file_record(record['path']) == record
            truth = pd.read_csv(inputs[f'{seed}:truth'], sep='\t')
            metadata = pd.read_csv(inputs[f'{seed}:metadata'], sep='\t')
            assert np.array_equal(truth.cell_id, metadata.cell_id)
            fit = pd.read_csv(folder/'fit/cells.tsv.gz', sep='\t')
            diff = pd.read_csv(folder/'fit/differences.tsv', sep='\t')
            pd.testing.assert_frame_equal(fit, pd.read_csv(inputs[f'{seed}:previous:cells.tsv.gz'], sep='\t'), check_exact=True)
            pd.testing.assert_frame_equal(diff, pd.read_csv(inputs[f'{seed}:previous:differences.tsv'], sep='\t'), check_exact=True)
            joined = fit.merge(truth[['cell_id', 'truth_clone']], on='cell_id', validate='one_to_one')
            s1 = joined[joined.specimen.eq('S1')]
            normal = s1.truth_clone.eq('normal')
            fpr = sum(s1.loc[normal, 'malignant_call'])/sum(normal)
            tpr = sum(s1.loc[~normal, 'malignant_call'])/sum(~normal)
            s2 = joined[joined.specimen.eq('S2') & joined.malignant_call]
            tumour = s1[s1.truth_clone.isin(['clone1', 'clone2']) & s1.malignant_call]
            groups = sorted(set(tumour.clone))
            table = [[int(sum((tumour.truth_clone == pop) & (tumour.clone == group))) for group in groups]
                for pop in ('clone1', 'clone2')]
            n = min(2, len(groups))
            accuracy = max(sum(table[i][i] for i in range(n)),
                sum(table[1-i][i] for i in range(n)))/len(tumour) if len(tumour) else None
            comparisons = diff[diff.specimen.eq('S1') & ~diff.clone_a.str.endswith('_c0') & ~diff.clone_b.str.endswith('_c0')]
            chrom = sorted(set(comparisons.chromosome))
            values = dict(normal_fpr=(fpr, fpr < .1), tumour_tpr=(tpr, tpr > .9),
                s2_one_group=(len(set(s2.clone)), len(set(s2.clone)) == 1),
                s1_two_groups=(len(groups), len(groups) == 2),
                assignment_accuracy=(accuracy, accuracy is not None and accuracy > .65),
                difference_chromosomes=(chrom, chrom == ['chr12']))
            old = json.loads(inputs[f'{seed}:previous:summary'].read_text())
            for name, (value, passed) in values.items():
                assert report['criteria'][name]['value'] == value == old['criteria'][name]['value']
                assert bool(report['criteria'][name]['passed']) == bool(passed) == old['criteria'][name]['passed']
            assert old['original_assertion_passed'] == all(v[1] for v in values.values())
            counts = Counter()
            for index, entry in enumerate(trace):
                assert index == entry['index'] and entry['complete']
                if entry['parent'] is not None:
                    assert 0 <= entry['parent'] < index
                saved = np.load(folder/entry['arrays'])
                rows = saved['rows']
                assert len(set(rows)) == len(rows)
                assert truth.iloc[rows].specimen.eq('S1').all()
                counts[entry['stage']] += 1
                if entry['stage'] == '_score_candidate':
                    cells, g = saved['cells'], saved['groups']
                    assert len(cells) == len(g) and len(set(cells)) == len(cells)
                    counts['accepted_candidates'] += int(entry['accepted'])
                    descendants = []
                    for item in trace[index+1:]:
                        parent = item['parent']
                        while parent is not None and parent != index:
                            parent = trace[parent]['parent']
                        if parent == index:
                            descendants.append(item)
                    agglomeration = [item for item in descendants if item['stage'] == 'agglomerate']
                    matched = [item for item in descendants if item['stage'] == 'matched_differences']
                    support = [item for item in descendants if item['stage'] == 'supported']
                    candidates.append(dict(seed=seed, index=index, source=entry['source'],
                        accepted=entry['accepted'], score=entry['score'], group_sizes=np.bincount(g).tolist(),
                        agglomerated_parts=[len(item['members']) for item in agglomeration],
                        matched_segment_counts=[[len(d['segments']) for d in item['differences']] for item in matched],
                        supported_counts=[len(item['kept']) for item in support]))
                elif entry['stage'] == 'event_fit' and entry['fitted']:
                    q = saved['q']; cells = saved['cells']; actual = truth.iloc[rows[cells]].truth_clone.to_numpy()
                    assert len(q) == len(cells) and np.isfinite(q).all() and np.all((q >= 0) & (q <= 1))
                    g = q > .5
                    tab = [[int(sum((actual == pop) & (g == group))) for group in (False, True)]
                        for pop in ('normal', 'clone1', 'clone2')]
                    num = sum(map(sum, tab[1:]))
                    accuracy_event = max(tab[1][0]+tab[2][1], tab[1][1]+tab[2][0])/num if num else None
                    events_out.append(dict(seed=seed, index=index, gain=entry['gain'],
                        cols=saved['cols'].tolist(), parameters=saved['parameters'].tolist(),
                        truth_table=tab, tumour_membership_accuracy=accuracy_event,
                        q_quantiles=np.quantile(q, [0,.25,.5,.75,1]).tolist()))
                elif entry['stage'] == 'matched_differences' and entry['caller'] == 'merge_unsupported':
                    labels, R = saved['labels'], saved['assignment_R']
                    assert np.array_equal(labels, R.argmax(axis=1))
                    assert R.shape == saved['full_R'].shape
                    assert np.isfinite(R).all() and (R >= 0).all()
                    # Summation budget from number of nonnegative probability terms.
                    budget = 8*R.shape[1]*np.finfo(R.dtype).eps
                    assert np.max(abs(R.sum(axis=1)-1)) <= budget
                    actual = truth.iloc[rows].truth_clone.to_numpy()
                    def contingency(labels):
                        return {pop: {str(k): int(sum((actual == pop) & (labels == k))) for k in sorted(set(labels))}
                            for pop in ('normal', 'clone1', 'clone2')}
                    memberships.append(dict(seed=seed, index=index, pairs=entry['pairs'],
                        proposal_assignment=contingency(labels),
                        full_assignment=contingency(saved['full_R'].argmax(axis=1)),
                        changed_assignments=int(sum(labels != saved['full_R'].argmax(axis=1))),
                        differences=entry['differences']))
                elif entry['stage'] == 'merge_unsupported':
                    before = np.load(folder/entry['input_arrays'])
                    assert np.array_equal(before['rows'], rows)
                    assert before['doses'].shape[0] >= saved['new_doses'].shape[0]
                    changed = before['doses'].shape[0] > saved['new_doses'].shape[0]
                    assert bool(entry['changed']) == changed
                    counts['changed_merges'] += int(changed)
            summaries.append(dict(seed=seed, original_criteria_pass=old['original_assertion_passed'],
                all_fit_tables_and_criteria_identical=True, captures=len(trace), stages=dict(counts)))
        write_json(dict(all_checks_passed=True, independent_criteria_checked=18, cases=summaries,
            candidates=candidates, merge_memberships=memberships, event_memberships=events_out,
            biological_specimens=0, production_changed=False, scientific_criteria_changed=False,
            interpretation='Known-noise observational trace, postfit truth descriptions only; no causal or calibration claim'),
            a.out/'report.json')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU independent criterion/identity/stage recount')


if __name__ == '__main__':
    main()
