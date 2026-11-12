"""Independent paired endpoints, membership matchings and simultaneous exact bounds."""
import argparse
from decimal import Decimal, localcontext
from functools import lru_cache
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta
from haplodose.provenance import start_run, finish_run, file_record

CONDITIONS = ('matched_nb', 'correlated_state', 'two_clones')
FAMILY = 168


@lru_cache(None)
def numerical_bound(k, n, direction, precision):
    if direction == 'upper' and k == n:
        return 1.
    if direction == 'lower' and k == 0:
        return 0.
    with localcontext() as context:
        context.prec = precision
        tail = Decimal(1)/Decimal(20*FAMILY)
        lo, hi = Decimal(0), Decimal(1)
        for _ in range(110):
            p = (lo+hi)/2
            odds, term, probability = p/(1-p), (1-p)**n, Decimal(0)
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


def bound(k, n, direction):
    alpha = .05/FAMILY
    scipy = (1. if k == n else float(beta.ppf(1-alpha, k+1, n-k))) if direction == 'upper' else (
             0. if k == 0 else float(beta.ppf(alpha, k, n-k+1)))
    first, second = (numerical_bound(k, n, direction, precision) for precision in (60, 80))
    assert first == second
    # 110 bisections and stable60/80-digit roots resolve far below float64 ULPs.
    assert abs(scipy-second) <= 128*np.finfo(float).eps
    return second


def verify_assignments(record):
    path = Path(record['assignments']['path'])
    assert file_record(path)['sha256'] == record['assignments']['sha256']
    with np.load(path) as z:
        truth, labels, eligible = z['truth'], z['raw_MAP'], z['eligible']
        assert len(truth) == len(labels) == 600
        assert len(set(eligible)) == len(eligible)
        table = [[int(sum((truth == k) & (labels == group))) for group in eligible] for k in (0, 1)]
        # Exhaustive matching is independent of the runner's Hungarian solver.
        if len(eligible) >= 2:
            best = max(table[0][i]+table[1][j] for i, j in itertools.permutations(range(len(eligible)), 2))
        elif len(eligible) == 1:
            best = max(table[0][0], table[1][0])
        else:
            best = 0
        accuracy = best/600
        values = dict(raw_candidate_count=len(eligible), multigroup=len(eligible) > 1,
            exactly_two=len(eligible) == 2, all_cell_one_to_one_accuracy=accuracy,
            correct_two=len(eligible) == 2 and accuracy > .65)
        for name, value in values.items():
            assert record[name] == value
    return values


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    old = a.task.parent/'dispersion_membership_20261109'
    inputs = dict(script=Path(__file__), design=a.task/'TRADEOFF_DESIGN.md',
        historical=old/'tradeoff_summary_v1/report.json')
    for c in CONDITIONS:
        for start in (0, 100):
            inputs[f'{c}_{start}_new'] = a.task/f'tradeoff_{c}_{start}_v1/panel/records.jsonl'
            inputs[f'{c}_{start}_old'] = old/f'tradeoff_MAP_{c}_{start}_v1/records.jsonl'
    rec = start_run(a.out, inputs, dict(family=FAMILY),
        dict(scope='Exposed development regression, biological n=0; no independent validation'), mode='offset_tradeoff_audit')
    try:
        summaries, pairs, cases = [], [], []
        for ci, c in enumerate(CONDITIONS):
            previous, current = [], []
            for start in (0, 100):
                folder = a.task/f'tradeoff_{c}_{start}_v1'
                for path in (folder/'run_manifest.json', folder/'panel/run_manifest.json'):
                    assert json.loads(path.read_text())['status'] == 'complete'
                prior = [json.loads(v) for v in inputs[f'{c}_{start}_old'].read_text().splitlines()]
                new = [json.loads(v) for v in inputs[f'{c}_{start}_new'].read_text().splitlines()]
                assert len(prior) == len(new) == 100
                for i, (x, y) in enumerate(zip(prior, new), start):
                    assert x['seed_key'] == y['seed_key'] == [20261109, 501, ci, i]
                    assert x['input_sha256'] == y['input_sha256']
                    assert x['condition'] == y['condition'] == c
                    verify_assignments(x)
                    values = verify_assignments(y)
                    cases.append(dict(condition=c, replicate=i, **values))
                previous.extend(prior)
                current.extend(new)
            metrics = ['multigroup'] if c != 'two_clones' else ['multigroup', 'exactly_two', 'correct_two']
            for metric in metrics:
                k = sum(bool(row[metric]) for row in current)
                direction = 'lower' if c == 'two_clones' else 'upper'
                value = bound(k, 200, direction)
                summaries.append(dict(condition=c, metric=metric, events=k, n=200, direction=direction,
                    family=FAMILY, simultaneous_bound=value,
                    error_criterion_met=value <= .05 if direction == 'upper' else None))
                pairs.append(dict(condition=c, metric=metric,
                    old_events=sum(bool(v[metric]) for v in previous), new_events=k,
                    lost=sum(bool(x[metric]) and not bool(y[metric]) for x, y in zip(previous, current)),
                    gained=sum(not bool(x[metric]) and bool(y[metric]) for x, y in zip(previous, current))))
        historical = json.loads(inputs['historical'].read_text())
        prior = historical['summary']+historical['retained_historical']
        assert historical['family'] == len(prior) == 163
        retained = []
        for row in prior:
            value = bound(row['events'], row['n'], row['direction'])
            retained.append(dict(row, family=FAMILY, simultaneous_bound=value,
                previous_bound=dict(family=163, value=row['simultaneous_bound']),
                error_criterion_met=value <= .05 if row['direction'] == 'upper' else None))
        assert len(cases) == 600 and len(summaries)+len(retained) == FAMILY
        (a.out/'report.json').write_text(json.dumps(dict(summary=summaries, retained_historical=retained,
            paired=pairs, cases=cases, family=FAMILY, paired_inputs=600, independent_assignment_matchings=1200,
            Decimal_precisions=[60, 80], all_checks_passed=True, biological_specimens=0,
            scientific_criteria_changed=False, scope='Known-reference binned development regression'), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU exhaustive assignment matching and Decimal60/80 exact binomial inversion')


if __name__ == '__main__':
    main()
