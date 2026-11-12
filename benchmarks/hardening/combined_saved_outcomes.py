"""Measure all saved combined outcomes without erasing failed preservation gates."""
import argparse
from itertools import permutations
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import file_record, start_run, finish_run


CASES = (('two_clone', 1), ('two_clone', 2), ('two_clone', 4), ('state', 3))


def assignments(table):
    """Exhaustive rectangular matching, independently of Hungarian scoring."""
    table = np.asarray(table, dtype=int)
    if not table.size:
        return 0
    if table.shape[0] > table.shape[1]:
        table = table.T
    return max(sum(int(table[row, col]) for row, col in enumerate(choice))
               for choice in permutations(range(table.shape[1]), table.shape[0]))


def score(folder, truth, condition):
    cells = pd.read_csv(folder/'cells.tsv.gz', sep='\t')
    cells = cells[cells.specimen.eq('S1')].merge(
        truth.loc[truth.specimen.eq('S1'), ['cell_id', 'truth_clone']],
        on='cell_id', how='outer', validate='one_to_one', indicator=True)
    assert len(cells) == 600 and cells._merge.eq('both').all()
    assert cells.malignant_call.dtype.kind == 'b'
    profiles = pd.read_csv(folder/'clones.tsv', sep='\t')
    assert profiles.normal_like.dtype.kind == 'b'
    names = sorted(profiles.loc[profiles.specimen.eq('S1') & ~profiles.normal_like, 'clone'])
    populations = sorted(cells.truth_clone.unique())
    table = np.array([[int((cells.truth_clone.eq(p) & cells.clone.eq(g)
                          & cells.malignant_call).sum()) for g in names] for p in populations])
    differences = pd.read_csv(folder/'differences.tsv', sep='\t')
    differences = differences[differences.specimen.eq('S1')
        & differences.clone_a.isin(names) & differences.clone_b.isin(names)]
    chromosomes = sorted(set(differences.chromosome))
    return dict(reported_groups=len(names), null_error=len(names) > 1,
        exact_chr12_recovery=len(names) == 2 and chromosomes == ['chr12'],
        difference_chromosomes=chromosomes, assignment_correct=assignments(table),
        assignment_n=600, tumour_calls=int(cells.malignant_call.sum()),
        assignment_table=table.tolist())


def input_control(before, after):
    count = 0
    for name in ('generation.npz',):
        with np.load(before/name) as x, np.load(after/name) as y:
            assert x.files == y.files
            for key in x.files:
                np.testing.assert_array_equal(x[key], y[key])
            count += 1
    with h5py.File(before/'input/counts.h5') as x, h5py.File(after/'input/counts.h5') as y:
        for key in ('X/data', 'X/indices', 'X/indptr', 'genes', 'cell_ids'):
            np.testing.assert_array_equal(x[key][:], y[key][:])
            count += 1
    for name in ('model_cells.tsv', 'truth_cells.tsv', 'positions.tsv'):
        pd.testing.assert_frame_equal(pd.read_csv(before/'input'/name, sep='\t'),
            pd.read_csv(after/'input'/name, sep='\t'), check_exact=True)
    a, b = before/'pooled_prevalence/folds', after/'pooled_prevalence/folds'
    names = sorted(p.name for p in a.glob('*.npz'))
    assert len(names) == 1 and names == sorted(p.name for p in b.glob('*.npz'))
    for name in names:
        with np.load(a/name) as x, np.load(b/name) as y:
            assert x.files == y.files
            for key in x.files:
                np.testing.assert_array_equal(x[key], y[key])
    assert count == 6
    return dict(generation_and_matrix_arrays_exact=6, metadata_tables_exact=3,
        candidate_fold_archives_exact=1, exact=True)


def output_control(before, after):
    eligible = lambda p: p.suffix == '.npy' or p.name.endswith(('.tsv', '.tsv.gz'))
    a, b = before/'production/fit', after/'production/fit'
    names = sorted(p.name for p in a.iterdir() if eligible(p))
    assert names == sorted(p.name for p in b.iterdir() if eligible(p))
    records = []
    for name in names:
        try:
            if name.endswith('.npy'):
                np.testing.assert_array_equal(np.load(a/name), np.load(b/name))
            else:
                pd.testing.assert_frame_equal(pd.read_csv(a/name, sep='\t'),
                    pd.read_csv(b/name, sep='\t'), check_exact=True)
            record = dict(file=name, exact=True)
        except AssertionError as error:
            record = dict(file=name, exact=False, discrepancy=str(error))
        records.append(record)
    assert len(records) == 16
    return dict(files=records, strict_original_preservation_passed=all(r['exact'] for r in records))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    inputs = dict(script=Path(__file__), design=args.task/'SAVED_COMBINED_READER_DESIGN.md',
        installed_wheel=args.task/'validation_23169639/installed_artifact.json')
    locations = []
    for condition, rep in CASES:
        seed = 202611098000+rep
        parent = args.task/f'combined_phase8_{condition}_{rep}_v1'
        before = args.task.parent/'pooled_prevalence_pilot_20261110'/(
            'pilot_'+condition+'_00_v1')/f'seed_{seed}'
        after = parent/'replay'/f'seed_{seed}'
        locations.append((condition, rep, parent, before, after))
        for key, path in dict(parent=parent/'run_manifest.json',
                failed_records=parent/'records_failed.json',
                report=parent/'report.json',
                rows=parent/'replay/results.json', old_rows=before.parent/'results.json',
                truth=after/'input/truth_cells.tsv').items():
            if key not in ('failed_records', 'report') or path.is_file():
                inputs[f'{condition}_{rep}_{key}'] = path
        for arm in ('production', 'pooled_prevalence'):
            for label, base in (('before', before), ('after', after)):
                folder = base/arm/'fit'
                inputs[f'{condition}_{rep}_{arm}_{label}_manifest'] = folder/'run_manifest.json'
                for name in ('cells.tsv.gz', 'clones.tsv', 'differences.tsv', 'log.tsv'):
                    inputs[f'{condition}_{rep}_{arm}_{label}_{name}'] = folder/name
        for label, base in (('before', before), ('after', after)):
            for name in ('generation.npz', 'input/counts.h5', 'input/model_cells.tsv',
                         'input/truth_cells.tsv', 'input/positions.tsv'):
                inputs[f'{condition}_{rep}_{label}_{name}'] = base/name
            for path in (base/'pooled_prevalence/folds').glob('*.npz'):
                inputs[f'{condition}_{rep}_{label}_fold_{path.name}'] = path
            for path in (base/'production/fit').iterdir():
                if path.suffix == '.npy' or path.name.endswith(('.tsv', '.tsv.gz')):
                    inputs[f'{condition}_{rep}_{label}_production_{path.name}'] = path
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, inputs, dict(cases=CASES, inference_rerun=False,
        original_strict_gates_preserved=True),
        dict(scope='Post-failure explanation of four exposed simulations; no revised promotion verdict'),
        mode='combined_saved_outcomes')
    try:
        # Independent exact assignment fixtures include both rectangular orientations.
        assert assignments([[324, 35, 1], [118, 100, 22]]) == 424
        assert assignments([[324, 36], [118, 122]]) == 446
        assert assignments([[360], [240]]) == 360
        assert assignments([[360, 240]]) == 360
        records = []
        for condition, rep, parent, before, after in locations:
            parent_manifest = json.loads((parent/'run_manifest.json').read_text())
            assert parent_manifest['status'] in ('complete', 'failed')
            if parent_manifest['status'] == 'complete':
                evidence = json.loads((parent/'report.json').read_text())
                support = evidence['returned_support_checks']
                previews = evidence['preview_records']
            else:
                evidence = json.loads((parent/'records_failed.json').read_text())
                support, previews = evidence['support'], evidence['previews']
            assert support and previews
            assert all(r['all_supported'] and r['shapes_aligned'] for r in support)
            for base in (before, after):
                for arm in ('production', 'pooled_prevalence'):
                    assert json.loads((base/arm/'fit/run_manifest.json').read_text())['status'] == 'complete'
            parity = input_control(before, after)
            preservation = output_control(before, after)
            current_rows = json.loads((parent/'replay/results.json').read_text())
            old_rows = [r for r in json.loads((before.parent/'results.json').read_text()) if r['replicate'] == rep]
            assert len(current_rows) == len(old_rows) == 2
            truth = pd.read_csv(after/'input/truth_cells.tsv', sep='\t')
            outcomes = []
            for label, base, rows in (('before', before, old_rows), ('after', after, current_rows)):
                for arm in ('production', 'pooled_prevalence'):
                    actual = score(base/arm/'fit', truth, condition)
                    official = next(r for r in rows if r['arm'] == arm)
                    assert all(actual[k] == official[k] for k in actual), (condition, rep, label, arm)
                    outcomes.append(dict(version=label, arm=arm, **actual))
            logs = pd.read_csv(after/'production/fit/log.tsv', sep='\t')
            dropped = logs.loc[logs.action.eq('drop_small_final')].dropna(axis=1, how='all')
            records.append(dict(condition=condition, replicate=rep,
                parent_status=parent_manifest['status'], exact_inputs=parity,
                production_preservation=preservation, outcomes=outcomes,
                returned_support_checks=len(support), preview_records=len(previews),
                all_recorded_returns_supported=True,
                final_support_removals=dropped.to_dict('records')))
        write_json(dict(cases=records, all_inputs_exact=True,
            independent_assignment_scores_match=True, original_strict_gates_unchanged=True,
            new_integration_pass_claimed=False, promotion=False, false_split_control_established=False,
            independent_evaluation=False, biological_specimens=0, frozen_results_changed=False),
            args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU saved identified outputs; exact independent assignment enumeration')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
