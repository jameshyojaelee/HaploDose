"""Reporting-only log-tail correction from immutable saved simulation outputs.

Validates all six coverage counts against the original interval table. It neither
reruns simulation/inference nor changes endpoints, hypotheses, intervals or data.
Underflowed p-values are null with explicit log probabilities, never literal zero.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import math
from pathlib import Path

from haplodose.io import new_output_dir, write_json
from haplodose.provenance import file_record, finish_run, start_run


def log_binomial_lower_tail(successes, total, probability):
    terms = [math.lgamma(total + 1) - math.lgamma(k + 1) - math.lgamma(total - k + 1)
             + k * math.log(probability) + (total - k) * math.log1p(-probability)
             for k in range(successes + 1)]
    largest = max(terms)
    return min(0., largest + math.log(math.fsum(math.exp(value - largest) for value in terms)))


def validate_counts(source, intervals):
    family = source['primary_undercoverage_family']
    if (family['endpoint'] != 'reported_mean_covers_zero' or family['size'] != 6
            or len(family['conditions']) != 3 or len(family['intervals']) != 2
            or family['family_alpha'] != .05):
        raise ValueError('Input is not the frozen six-member primary undercoverage family')
    expected = {(condition, interval) for condition in family['conditions'] for interval in family['intervals']}
    grouped = {}
    with intervals.open(newline='') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            key = row['condition'], row['interval']
            if key not in expected:
                raise ValueError(f'Unexpected interval-table family member: {key}')
            truth = row['reported_mean_covers_zero']
            if truth not in ('True', 'False', '1', '0'):
                raise ValueError(f'Unrecognized coverage indicator: {truth}')
            group = grouped.setdefault(key, {})
            if row['replicate'] in group:
                raise ValueError(f'Duplicate replicate in primary family member: {key}')
            group[row['replicate']] = truth in ('True', '1')
    if set(grouped) != expected or len(source['summary']) != len(expected):
        raise ValueError('Incomplete primary family')
    seen = set()
    for row in source['summary']:
        key = row['condition'], row['interval']
        if key in seen:
            raise ValueError('Duplicate summary family member')
        seen.add(key)
        recorded = row['reported_mean_covers_zero']
        if (recorded['n_simulations'] != len(grouped[key]) or recorded['successes'] != sum(grouped[key].values())
                or recorded['n_simulations'] != source['replicates']):
            raise ValueError(f'Summary counts disagree with saved intervals: {key}')
    return grouped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task-root', type=Path, required=True)
    parser.add_argument('--source-summary', type=Path, required=True)
    parser.add_argument('--source-intervals', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root, out = args.task_root.resolve(strict=True), args.out.resolve()
    relative = out.relative_to(root)
    if not relative.parts or relative.parts[0] == 'baseline' or out.exists():
        raise ValueError('Output must be a fresh task_root descendant outside baseline')
    out = new_output_dir(out)
    original_manifest = args.source_summary.parent / 'run_manifest.json'
    record = start_run(out, dict(probe_source=__file__, original_summary=args.source_summary,
        original_intervals=args.source_intervals, original_manifest=original_manifest),
        dict(operation='reporting_only_log_binomial_tails', family_size=6, family_alpha=.05),
        dict(scope='same saved simulation counts and predeclared primary family',
             biological_replication='none'), mode='uncertainty_reporting_revision')
    try:
        source = json.loads(args.source_summary.read_text())
        original = json.loads(original_manifest.read_text())
        if original['status'] != 'complete':
            raise ValueError('Original simulation manifest is not complete')
        validate_counts(source, args.source_intervals)
        corrected = copy.deepcopy(source)
        family = source['primary_undercoverage_family']
        target, alpha, size = family['nominal_coverage_target'], family['family_alpha'], family['size']
        logs = [log_binomial_lower_tail(row['reported_mean_covers_zero']['successes'],
                row['reported_mean_covers_zero']['n_simulations'], target) for row in source['summary']]
        order = sorted(range(size), key=lambda k: logs[k])
        adjusted, previous = [0.] * size, -math.inf
        for rank, index in enumerate(order):
            previous = max(previous, logs[index] + math.log(size - rank))
            adjusted[index] = min(0., previous)
        smallest_log = math.log(math.ulp(0.))
        table = []
        for row, log_p, log_holm in zip(corrected['summary'], logs, adjusted):
            primary = row['reported_mean_covers_zero']
            exact_underflow, holm_underflow = log_p < smallest_log, log_holm < smallest_log
            primary.update(original_reported_p_exact=primary['p_undercoverage_exact'],
                original_reported_p_holm=primary['p_undercoverage_holm'],
                p_undercoverage_exact=None if exact_underflow else math.exp(log_p),
                p_undercoverage_holm=None if holm_underflow else math.exp(log_holm),
                log_p_undercoverage_exact=log_p, log_p_undercoverage_holm=log_holm,
                log10_p_undercoverage_exact=log_p / math.log(10),
                log10_p_undercoverage_holm=log_holm / math.log(10),
                exact_p_underflow=exact_underflow, holm_p_underflow=holm_underflow,
                p_value_representation='null means below float64 range; use log probabilities',
                falsified_undercoverage=log_holm < math.log(alpha))
            table.append(dict(condition=row['condition'], interval=row['interval'],
                successes=primary['successes'], n_simulations=primary['n_simulations'],
                log10_p_exact=primary['log10_p_undercoverage_exact'], log10_p_holm=primary['log10_p_undercoverage_holm'],
                simultaneous_cp_ci_low=primary['simultaneous_cp_ci_low'],
                simultaneous_cp_ci_high=primary['simultaneous_cp_ci_high'],
                falsified_undercoverage=primary['falsified_undercoverage']))
        corrected['reporting_revision'] = dict(operation='log-tail reporting correction only',
            count_validation='all six counts matched immutable saved interval records',
            source_summary=file_record(args.source_summary), source_intervals=file_record(args.source_intervals),
            source_manifest=file_record(original_manifest), original_output_id=original['output']['id'],
            tail_method='independent math.lgamma log-PMF summation with log-sum-exp; Holm in log space',
            unchanged=['simulation counts', 'nominal coverage target', 'primary hypothesis family',
                       'simultaneous CP intervals', 'all descriptive endpoints', 'biological replication status'])
        write_json(corrected, out / 'corrected_summary.json')
        with (out / 'primary_logtails.tsv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(table[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(table)
        finish_run(out, record, backend='reporting_only_no_inference')
        print(json.dumps(dict(output=str(out / 'corrected_summary.json'),
                              falsified_primary_members=sum(row['falsified_undercoverage'] for row in table))))
    except Exception as error:
        finish_run(out, record, error=error)
        raise


if __name__ == '__main__':
    main()
