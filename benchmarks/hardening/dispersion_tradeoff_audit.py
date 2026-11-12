"""Independent endpoint matching and high-precision binomial-bound inversion."""
import argparse
from decimal import Decimal, localcontext
from functools import lru_cache
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.provenance import start_run, finish_run, file_record


@lru_cache(None)
def independent_bound(k, n, direction, precision):
    if direction == 'upper' and k == n:
        return 1.
    if direction == 'lower' and k == 0:
        return 0.
    with localcontext() as context:
        context.prec = precision
        tail = Decimal(1)/Decimal(20*163)
        lo, hi = Decimal(0), Decimal(1)
        for _ in range(110):
            p = (lo+hi)/2
            odds = p/(1-p)
            term = (1-p)**n
            probability = Decimal(0)
            for j in range(n+1):
                if (j <= k if direction == 'upper' else j >= k):
                    probability += term
                if j < n:
                    term = term*(n-j)/(j+1)*odds
            if (probability > tail if direction == 'upper' else probability < tail):
                lo = p
            else:
                hi = p
        return float((lo+hi)/2)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    folder = a.task/'tradeoff_summary_v1'
    run = start_run(a.out, dict(script=Path(__file__), design=a.task/'TRADEOFF_DESIGN.md',
        summary=folder/'report.json', events=folder/'events.tsv'), {},
        dict(scope='Independent endpoints and all163 high-precision bounds; biological n=0'),
        mode='dispersion_tradeoff_audit')
    try:
        report = json.loads((folder/'report.json').read_text())
        events = pd.read_csv(folder/'events.tsv', sep='\t')
        assert report['family'] == 163 and len(events) == 1200
        bounds = report['summary']+report['retained_historical']
        assert len(bounds) == 163
        for row in bounds:
            first = independent_bound(row['events'], row['n'], row['direction'], 60)
            second = independent_bound(row['events'], row['n'], row['direction'], 80)
            assert first == second
            # Precision and root convergence define this representable-number budget.
            budget = 128*np.finfo(float).eps*max(1., abs(second))
            assert abs(row['simultaneous_bound']-second) <= budget
            if row['direction'] == 'upper':
                assert bool(row['error_criterion_met']) == (second <= .05)
        previous_path = a.task.parent/'joint_fresh_null_20261108/summary_v1/summary.json'
        previous = json.loads(previous_path.read_text())
        original = previous['summary']+previous['retained_historical']
        for old, new in zip(original, report['retained_historical']):
            assert (old['events'], old['n'], old['direction']) == (new['events'], new['n'], new['direction'])
            assert new['previous_bound']['value'] == old['simultaneous_bound']
        paired = {}
        for arm in ('MAP', 'soft'):
            for condition in ('matched_nb', 'correlated_state', 'two_clones'):
                cases = []
                for start in (0, 100):
                    directory = a.task/f'tradeoff_{arm}_{condition}_{start}_v1'
                    manifest = json.loads((directory/'run_manifest.json').read_text())
                    assert manifest['status'] == 'complete'
                    records = [json.loads(line) for line in (directory/'records.jsonl').read_text().splitlines()]
                    assert len(records) == 100
                    for record in records:
                        assert file_record(record['assignments']['path'])['sha256'] == record['assignments']['sha256']
                        with np.load(record['assignments']['path']) as z:
                            truth, labels, eligible = z['truth'], z['raw_MAP'], z['eligible']
                            assert len(truth) == len(labels) == 600
                            assert np.array_equal(labels, z['raw_R'].argmax(axis=1))
                            sizes = {int(k): int(np.sum(labels == k)) for k in np.unique(labels)}
                            expected_eligible = [k for k in range(1, z['raw_R'].shape[1]) if sizes.get(k, 0) >= 30]
                            assert expected_eligible == eligible.tolist()
                            best = max(sum(int(np.sum((truth == t) & (labels == k))) for t, k in enumerate(pair))
                                for pair in itertools.permutations(eligible.tolist()+[-2, -3], 2))
                            accuracy = best/600
                            assert accuracy == record['all_cell_one_to_one_accuracy']
                            assert bool(record['multigroup']) == (len(eligible) > 1)
                            assert bool(record['exactly_two']) == (len(eligible) == 2)
                            assert bool(record['correct_two']) == (len(eligible) == 2 and accuracy > .65)
                        cases.append(record)
                assert len(cases) == 200
                paired[(arm, condition)] = cases
                table = events[(events.arm == arm) & (events.condition == condition)]
                assert len(table) == 200
                for metric in ('multigroup', 'exactly_two', 'correct_two'):
                    assert sum(bool(v[metric]) for v in cases) == int(table[metric].sum())
                for row in report['summary']:
                    if row['arm'] == arm and row['condition'] == condition:
                        assert row['events'] == sum(bool(v[row['metric']]) for v in cases)
        for condition in ('matched_nb', 'correlated_state', 'two_clones'):
            for x, y in zip(paired[('MAP', condition)], paired[('soft', condition)]):
                assert x['input_sha256'] == y['input_sha256']
        for row in report['paired']:
            x, y = paired[('MAP', row['condition'])], paired[('soft', row['condition'])]
            metric = row['metric']
            assert row['lost'] == sum(bool(a[metric]) and not bool(b[metric]) for a, b in zip(x, y))
            assert row['gained'] == sum(not bool(a[metric]) and bool(b[metric]) for a, b in zip(x, y))
        (a.out/'report.json').write_text(json.dumps(dict(fits=1200, independently_matched_assignments=1200,
            binomial_bounds=163, decimal_precision=[60, 80], historical_endpoints_preserved=153,
            biological_specimens=0, all_checks_passed=True, production_changed=False,
            scientific_criteria_changed=False), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, run, error=error)
        raise
    finish_run(a.out, run, backend='CPU exhaustive label matching, scalar endpoint recount, Decimal60/80 polynomial inversion')


if __name__ == '__main__':
    main()
