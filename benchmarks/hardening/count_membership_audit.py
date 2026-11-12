"""Independently recount count-oracle and saved-partition endpoints."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    source = args.task/'reference_v2'
    inputs = dict(script=Path(__file__), manifest=source/'run_manifest.json')
    manifest = json.loads(inputs['manifest'].read_text())
    assert manifest['status'] == 'complete'
    for group in ('input','output'):
        for record in manifest[group]['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
    for name in ('cases.tsv','posterior.npz','report.json','quadrature.json','sufficiency_checks.json'):
        inputs[name] = source/name
    run = start_run(args.out,inputs,dict(cases=200),
                    dict(scope='Independent endpoint recount; no new method or biological evidence'),
                    mode='count_membership_audit')
    try:
        frame = pd.read_csv(inputs['cases.tsv'],sep='\t')
        assert frame.replicate.tolist() == list(range(200))
        with np.load(inputs['posterior.npz']) as data:
            probability = data['known_depth_count_probability']
            convergence = data['max_log_odds_changes']
        assert probability.shape == (200,600) and convergence.shape == (200,4)
        assert np.isfinite(probability).all() and ((probability>=0)&(probability<=1)).all()
        assert np.isfinite(convergence).all() and (convergence>=0).all()
        checked = 0
        for i, row in frame.iterrows():
            truth = np.repeat([0,1],300)
            assert row.known_depth_count_correct == int(np.sum((probability[i]>.5)==truth))
            balanced = np.zeros(600,int)
            # Sorting probabilities is monotonic in odds; ties could lose
            # ordering information, so require no tie at the boundary.
            order = np.argsort(probability[i],kind='stable')
            assert probability[i,order[299]] < probability[i,order[300]]
            balanced[order[300:]] = 1
            assert row.known_depth_balanced_correct == int(np.sum(balanced==truth))
            for arm in ('joint_full','joint_proposal'):
                base = args.task.parents[1]/'dev/joint_fresh_null_20261108/chunks_v1'/arm/f'eval_two_{i//50*50}'
                path = base/f'two_clones_120_{i:03d}.npz'
                assert file_record(path)['sha256'] == row[arm+'_assignment_sha256']
                with np.load(path) as data:
                    assert np.array_equal(data['simulated_CNA_labels'],truth)
                    labels = data['raw_MAP']
                    groups = data['eligible_tumour_groups']
                    table = [[int(np.sum((truth==t)&(labels==g))) for g in groups] + [0,0] for t in (0,1)]
                    best = max(table[0][a]+table[1][b] for a,b in itertools.permutations(range(len(groups)+2),2))
                    many = sum(max(table[0][g],table[1][g]) for g in range(len(groups)))
                    assert best == row[arm+'_correct']
                    assert many == row[arm+'_oracle_coarsening_correct']
                    amplitude = data['latent_state']+(.6*np.log(1.5))*truth
                    assert row.noiseless_amplitude_correct == int(np.sum((amplitude>.3*np.log(1.5))==truth))
                    checked += 1
        report = json.loads(inputs['report.json'].read_text())
        for name,value in report['mean_accuracy'].items():
            assert float(frame[name].mean()/600) == value
        assert np.array_equal(convergence.max(axis=0),report['max_log_odds_changes_by_order'])
        assert frame.count_decisions_changed_128_to_256.sum() == report['decisions_changed_128_to_256']
        quadrature = json.loads(inputs['quadrature.json'].read_text())
        assert {r['replicate'] for r in quadrature} == {0,99,199}
        assert all(np.isfinite(r['odds_difference']) and all(0<=v<1 for v in r['prior_tail_relative_bound']) for r in quadrature)
        assert max(abs(r['odds_difference']) for r in quadrature) == report['max_independent_log_odds_discrepancy']
        (args.out/'report.json').write_text(json.dumps(dict(cases_recounted=200,
            independently_enumerated_assignment_matchings=checked, conditional_oracle_only=True,
            biological_validation=False, scientific_criteria_changed=False),indent=2))
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error)
        raise


if __name__ == '__main__':
    main()
