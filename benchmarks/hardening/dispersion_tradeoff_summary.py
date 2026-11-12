"""Collect every paired dispersion fit and retain the expanded confidence family."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta
from haplodose.provenance import start_run, finish_run, file_record

ARMS = ('MAP', 'soft')
CONDITIONS = ('matched_nb', 'correlated_state', 'two_clones')
FAMILY = 163


def bound(k, n, direction):
    alpha = .05/FAMILY
    if direction == 'upper':
        return 1. if k == n else float(beta.ppf(1-alpha, k+1, n-k))
    return 0. if k == 0 else float(beta.ppf(alpha, k, n-k+1))


def collect(task):
    result, inputs = {}, {}
    for arm in ARMS:
        for condition in CONDITIONS:
            rows = []
            for start in (0, 100):
                folder = task/f'tradeoff_{arm}_{condition}_{start}_v1'
                manifest = json.loads((folder/'run_manifest.json').read_text())
                if manifest['status'] != 'complete':
                    raise ValueError('Incomplete scientific chunk')
                inputs[f'{arm}_{condition}_{start}_manifest'] = folder/'run_manifest.json'
                inputs[f'{arm}_{condition}_{start}_records'] = folder/'records.jsonl'
                block = [json.loads(line) for line in (folder/'records.jsonl').read_text().splitlines()]
                if len(block) != 100:
                    raise ValueError('Missing or duplicate fits')
                for i, row in enumerate(block, start):
                    if row['seed_key'] != [20261109, 501, CONDITIONS.index(condition), i]:
                        raise ValueError('Seed order changed')
                    if row['arm'] != arm or row['condition'] != condition:
                        raise ValueError('Wrong scientific arm')
                    assert file_record(row['assignments']['path'])['sha256'] == row['assignments']['sha256']
                rows.extend(block)
            result[(arm, condition)] = rows
    for condition in CONDITIONS:
        for first, second in zip(result[('MAP', condition)], result[('soft', condition)]):
            if first['input_sha256'] != second['input_sha256']:
                raise ValueError('Unpaired observations')
    return result, inputs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    rows, inputs = collect(a.task)
    historical_path = a.task.parent/'joint_fresh_null_20261108/summary_v1/summary.json'
    inputs.update(script=Path(__file__), design=a.task/'TRADEOFF_DESIGN.md', historical=historical_path)
    rec = start_run(a.out, inputs, dict(confidence_family=FAMILY),
        dict(scope='Conditional bin workflow, biological n=0; no default promotion'), mode='dispersion_tradeoff_summary')
    try:
        summaries, paired, event_rows = [], [], []
        for arm in ARMS:
            for condition in CONDITIONS:
                cases = rows[(arm, condition)]
                metrics = ['multigroup'] if condition != 'two_clones' else ['multigroup', 'exactly_two', 'correct_two']
                for metric in metrics:
                    events = sum(bool(v[metric]) for v in cases)
                    direction = 'upper' if condition != 'two_clones' else 'lower'
                    value = bound(events, 200, direction)
                    summaries.append(dict(arm=arm, condition=condition, metric=metric, events=events, n=200,
                        direction=direction, simultaneous_bound=value, family=FAMILY,
                        error_criterion_met=bool(value <= .05) if direction == 'upper' else None))
                for v in cases:
                    event_rows.append({key: v[key] for key in ('arm', 'condition', 'replicate', 'multigroup',
                        'exactly_two', 'correct_two', 'all_cell_one_to_one_accuracy', 'raw_candidate_count')})
        for condition in CONDITIONS:
            metrics = ['multigroup'] if condition != 'two_clones' else ['multigroup', 'exactly_two', 'correct_two']
            for metric in metrics:
                x, y = rows[('MAP', condition)], rows[('soft', condition)]
                paired.append(dict(condition=condition, metric=metric,
                    lost=sum(bool(a[metric]) and not bool(b[metric]) for a, b in zip(x, y)),
                    gained=sum(not bool(a[metric]) and bool(b[metric]) for a, b in zip(x, y))))
        old = json.loads(historical_path.read_text())
        previous = old['summary']+old['retained_historical']
        assert old['family'] == 153 and len(previous) == 153 and len(summaries) == 10
        retained = []
        for row in previous:
            value = bound(row['events'], row['n'], row['direction'])
            updated = dict(row)
            updated.update(family=FAMILY, simultaneous_bound=value,
                previous_bound=dict(family=153, value=row['simultaneous_bound']),
                error_criterion_met=bool(value <= .05) if row['direction'] == 'upper' else None)
            retained.append(updated)
        pd.DataFrame(summaries).to_csv(a.out/'summary.tsv', sep='\t', index=False)
        pd.DataFrame(event_rows).to_csv(a.out/'events.tsv', sep='\t', index=False)
        (a.out/'report.json').write_text(json.dumps(dict(summary=summaries, retained_historical=retained,
            paired=paired, family=FAMILY, total_fits=1200, paired_observations=600, biological_specimens=0,
            production_changed=False, scientific_criteria_changed=False,
            scope='Known-reference bin-only workflow; no biological or public raw-input error guarantee'), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU complete paired collection and exact binomial family bounds')


if __name__ == '__main__':
    main()
