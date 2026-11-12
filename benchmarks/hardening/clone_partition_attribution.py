"""Descriptive attribution of saved clone fragmentation and membership mixing."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from haplodose.provenance import start_run, finish_run, file_record


def verify(record):
    assert file_record(record['path'])['sha256'] == record['sha256'], record['path']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    source = args.task.parents[1] / 'dev/joint_fresh_null_20261108'
    inputs = dict(script=Path(__file__), design=args.task / 'DESIGN.md',
                  prior_audit=source / 'audit_v1/run_manifest.json')
    assert json.loads(inputs['prior_audit'].read_text())['status'] == 'complete'
    paths = []
    for arm, role, start in itertools.product(('joint_proposal', 'joint_full'),
                                             ('eval_state', 'eval_two'), (0, 50, 100, 150)):
        directory = source / 'chunks_v1' / arm / f'{role}_{start}'
        manifest = json.loads((directory / 'run_manifest.json').read_text())
        assert manifest['status'] == 'complete'
        verify(manifest['output']['files']['records.jsonl'])
        key = f'{arm}:{role}:{start}'
        inputs[key + ':manifest'] = directory / 'run_manifest.json'
        inputs[key + ':records'] = directory / 'records.jsonl'
        paths.append((arm, role, start, inputs[key + ':records']))
    run = start_run(args.out, inputs, dict(cases=800, refits=0),
                    dict(scope='Descriptive simulated latent attribution only'),
                    mode='clone_partition_attribution')
    try:
        records = []
        for arm, role, start, path in paths:
            seen = []
            with path.open() as stream:
                for line in stream:
                    row = json.loads(line)
                    seen.append(row['replicate'])
                    assert row['seed_key'] == [20261107, 120, 1 if role == 'eval_state' else 2, row['replicate']]
                    verify(row['assignments'])
                    with np.load(row['assignments']['path']) as data:
                        truth = data['simulated_CNA_labels']
                        labels = data['raw_MAP']
                        eligible = data['eligible_tumour_groups']
                        state = data['latent_state']
                    assert truth.shape == labels.shape == (600,)
                    if state.ndim == 2:
                        # Generator stores a cell-by-bin regional effect; use
                        # all coordinates, avoiding a selected genomic region.
                        assert state.shape[0] == len(truth)
                    else:
                        assert state.shape == truth.shape
                        state = state[:, None]
                    assert np.isfinite(state).all()
                    actual = np.unique(truth)
                    table = np.array([[np.sum((truth == t) & (labels == g)) for g in eligible]
                                      for t in actual], dtype=int).reshape(len(actual), len(eligible))
                    covered = int(table.sum())
                    one = int(table[linear_sum_assignment(-table)].sum()) if table.size else 0
                    many = int(table.max(axis=0).sum()) if table.size else 0
                    # Independently enumerate injective matching for two true
                    # clones, including an unassigned dummy output group.
                    assert len(actual) in (1, 2)
                    padded = np.pad(table, ((0, 0), (0, len(actual))))
                    exhaustive = max(sum(padded[i, g] for i, g in enumerate(groups))
                                     for groups in itertools.permutations(range(padded.shape[1]), len(actual)))
                    assert one == exhaustive and 0 <= one <= many <= covered <= 600
                    entropy = 0.
                    for column in table.T:
                        if not column.sum():
                            continue
                        probability = column[column > 0] / column.sum()
                        entropy -= column.sum() * float(np.sum(probability * np.log2(probability))) / 600
                    residual = np.empty_like(state, dtype=float)
                    fitted = np.zeros_like(state, dtype=float)
                    for t in actual:
                        mask = truth == t
                        residual[mask] = state[mask] - state[mask].mean(axis=0)
                        for g in np.unique(labels[mask]):
                            selected = mask & (labels == g)
                            fitted[selected] = residual[selected].mean(axis=0)
                    denominator = float(np.sum(residual ** 2))
                    explained = float(np.sum(fitted ** 2)) / denominator if denominator else np.nan
                    records.append(dict(arm=arm, role=role, replicate=row['replicate'],
                        groups=len(eligible), coverage=covered/600, one_to_one=one/600,
                        oracle_many_to_one=many/600, fragmentation_gap=(many-one)/600,
                        truth_entropy_covered_weighted=entropy,
                        largest_eligible_fraction=int(table.sum(axis=0).max())/600 if table.size else 0.,
                        state_variance_explained_within_truth=explained,
                        assignment_sha256=row['assignments']['sha256']))
            assert seen == list(range(start, start+50))
        frame = pd.DataFrame(records)
        assert len(frame) == 800
        frame.to_csv(args.out / 'cases.tsv', sep='\t', index=False)
        metrics = ['groups', 'coverage', 'one_to_one', 'oracle_many_to_one',
                   'fragmentation_gap', 'truth_entropy_covered_weighted',
                   'state_variance_explained_within_truth']
        frame.groupby(['arm', 'role'])[metrics].agg(['mean', 'median', 'min', 'max', 'count']).to_csv(
            args.out / 'descriptive_summary.tsv', sep='\t')
        (args.out / 'report.json').write_text(json.dumps(dict(cases=800, refits=0,
            independent_assignment_matching_checks=800, scientific_gate_changed=False,
            biological_validation=False), indent=2))
        for record in run['input']['files'].values():
            verify(record)
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
