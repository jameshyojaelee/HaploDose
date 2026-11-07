"""Summarize frozen clone diagnostics without refitting or adding replicates.

Report both saved separate-family bounds and a single one-sided 95% family of nine bounds:
six upper bounds for null extra-group rates plus three lower bounds for same-nuisance alternative
multigroup detection. Detection power does not mean correct CNA membership. Run on a compute node.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from importlib import metadata
import json
from pathlib import Path
import platform
import sys

import numpy as np
from scipy import stats


METHODS = ('raw', 'selected_filter', 'full_replay')


def file_record(path):
    path = Path(path).resolve()
    before = path.stat()
    content = path.read_bytes()
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError(f'Input changed while reading: {path}')
    return dict(path=str(path), bytes=after.st_size, sha256=hashlib.sha256(content).hexdigest())


def distribution(values):
    values = np.asarray(values, float)
    if not len(values) or not np.isfinite(values).all():
        raise ValueError('Need finite nonempty saved ARI values')
    return dict(n=len(values), mean=float(values.mean()), median=float(np.median(values)),
                quartile1=float(np.quantile(values, .25)), quartile3=float(np.quantile(values, .75)),
                minimum=float(values.min()), maximum=float(values.max()))


def summarize(directory, condition):
    directory = Path(directory)
    files = [directory / name for name in ('provenance.json', 'inputs.jsonl', 'replicates.jsonl', 'summary.json')]
    manifest = [file_record(path) for path in files]
    provenance = json.loads(files[0].read_text())
    inputs = [json.loads(line) for line in files[1].read_text().splitlines()]
    rows = [json.loads(line) for line in files[2].read_text().splitlines()]
    saved = json.loads(files[3].read_text())
    n = provenance['design']['replicates']
    if len(rows) != n or len(inputs) != n or sorted(row['replicate'] for row in rows) != list(range(n)):
        raise ValueError(f'{condition}: incomplete or duplicated replicate records')
    if any(row['condition'] != condition for row in rows):
        raise ValueError(f'{condition}: mismatched condition label')
    by_input = {row['replicate']: row for row in inputs}
    if len(by_input) != n or any(row['input_sha256'] != by_input[row['replicate']]['input_sha256'] for row in rows):
        raise ValueError(f'{condition}: input hashes do not match replicate records')
    if condition == 'two_clones':
        if (saved.get('generator_version') != 'same_nuisance_v2'
                or any(row.get('generator_version') != 'same_nuisance_v2' for row in rows)
                or set(saved.get('nuisance_model', {}).get('applies_to', [])) != {'correlated_state', 'two_clones'}):
            raise ValueError('Primary alternative must be the corrected same_nuisance_v2 generator')
    endpoints = []
    for method in METHODS:
        events = [bool(row['full_replay_detected']) if method == 'full_replay'
                  else row[method + '_clones'] > 1 for row in rows]
        successes = sum(events)
        original = next(row for row in saved['summary'] if row['method'] == method)
        if original['detections'] != successes or original['n'] != n:
            raise ValueError(f'{condition}/{method}: saved endpoint disagrees with replicate records')
        numerical_check = None
        if condition == 'two_clones':
            bound = 0. if successes == 0 else float(stats.beta.ppf(.05 / 9, successes, n - successes + 1))
            direction = 'lower'
        else:
            bound = 1. if successes == n else float(stats.beta.ppf(1 - .05 / 9, successes + 1, n - successes))
            direction = 'upper'
            if successes == 0:
                closed_form = float(-np.expm1(np.log(.05 / 9) / n))
                beta_decision, closed_decision = bound <= .05, closed_form <= .05
                numerical_check = dict(primary='closed_form', closed_form=closed_form, beta_ppf=bound,
                    absolute_discrepancy=abs(bound - closed_form), threshold=.05,
                    closed_form_margin_to_threshold=closed_form - .05, beta_ppf_margin_to_threshold=bound - .05,
                    threshold_decision=('unresolved' if beta_decision != closed_decision else
                                        'at_or_below_5pct' if closed_decision else 'above_5pct'),
                    rule='Compare both 5% decisions; discrepancies crossing the threshold remain unresolved')
                bound = closed_form
        endpoints.append(dict(method=method, events=successes, n=n, rate=successes / n,
                              joint_family9_bound_direction=direction, joint_family9_bound=bound,
                              zero_event_numerical_check=numerical_check,
                              saved_separate_family_summary=original))
    memberships = {}
    for method in ('raw', 'selected_filter'):
        counts = [row[method + '_clones'] for row in rows]
        memberships[method] = dict(ARI_against_simulated_CNA_labels=distribution([row[method + '_membership_ari'] for row in rows]),
                                   counted_tumour_group_distribution=dict(sorted(Counter(counts).items())),
                                   exactly_two_counted_groups=sum(value == 2 for value in counts),
                                   more_than_two_counted_groups=sum(value > 2 for value in counts),
                                   fewer_than_two_counted_groups=sum(value < 2 for value in counts))
    paired = dict(raw_selected_count_changes=sum(row['raw_clones'] != row['selected_filter_clones'] for row in rows),
                  raw_replay_event_changes=sum((row['raw_clones'] > 1) != row['full_replay_detected'] for row in rows),
                  selected_minus_raw_ARI=distribution([row['selected_filter_membership_ari'] - row['raw_membership_ari']
                                                       for row in rows]))
    # Verify the four input files still have their original content after summarization.
    if [file_record(path) for path in files] != manifest:
        raise RuntimeError(f'{condition}: source records changed during summary')
    return dict(condition=condition, endpoints=endpoints, membership=memberships, paired=paired,
                input_manifest=manifest, generator_version=saved.get('generator_version', 'v4_original'),
                membership_scope='all-cell MAP labels including reference-like bucket; no explicit cell abstention')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matched', type=Path, required=True)
    parser.add_argument('--state', type=Path, required=True)
    parser.add_argument('--power', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    results = [summarize(path, condition) for path, condition in
               ((args.matched, 'matched_nb'), (args.state, 'correlated_state'), (args.power, 'two_clones'))]
    output = dict(results=results, confidence=dict(joint_family_size=9, one_sided_family_alpha=.05,
                      per_bound_alpha=.05 / 9, method='Bonferroni exact binomial',
                      endpoints='six null-rate upper bounds and three alternative detection lower bounds',
                      saved_families='separate null-family6 and power-family3 summaries retained unchanged',
                      zero_event_upper_at_n100=float(-np.expm1(np.log(.05 / 9) / 100))),
                  limitations=['simulation specimens, not independent patients',
                               'known reference; binned expression only, no gene veto, alleles or spatial graph',
                               'multigroup detection does not establish correct CNA cell membership',
                               'full replay emits a global gate, not new membership or cell abstention',
                               'constant one-clone truth makes ARI degenerate; it is not cell accuracy',
                               'no fixed scientific minimum power was specified; no validated-clone release claim'],
                  validation_status='diagnostic_only', provenance=dict(source=file_record(__file__), command=sys.argv,
                      executable=sys.executable, python=sys.version, host=platform.node(),
                      versions={name: metadata.version(name) for name in ('numpy', 'scipy')}))
    (args.out / 'summary.json').write_text(json.dumps(output, indent=2) + '\n')


if __name__ == '__main__':
    main()
