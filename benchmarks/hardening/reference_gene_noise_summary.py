"""Collect paired pilot error and power at the frozen scientific operating point."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import beta
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'PILOT_DESIGN.md')
    for condition in ('matched', 'state', 'two_clone'):
        for start in (0, 5):
            inputs[f'{condition}{start}'] = args.task/f'pilot_{condition}_{start:02d}_v1/results.json'
    rec = start_run(args.out, inputs, dict(error_family='2 nulls x2 looks', alpha_per_look=.0125),
        dict(scope='Small fresh conditional workflow pilot, not biological validation'), mode='reference_gene_noise_summary')
    try:
        results = []
        for name, path in inputs.items():
            if name not in ('script', 'design'):
                results.extend(json.loads(path.read_text()))
        assert len(results) == 60
        summaries = []
        for condition in ('matched', 'state', 'two_clone'):
            for arm in ('production', 'candidate'):
                selected = [r for r in results if r['condition'] == condition and r['arm'] == arm]
                assert len(selected) == 10 and sorted(r['replicate'] for r in selected) == list(range(10))
                n = len(selected)
                errors = sum(r['null_error'] for r in selected)
                recovery = sum(r['exact_chr12_recovery'] for r in selected)
                null = condition != 'two_clone'
                probability = Fraction(sum(math.comb(n, j)*19**(n-j) for j in range(errors, n+1)), 20**n)
                upper = 1. if errors == n else float(beta.ppf(1-.0125, errors+1, n-errors))
                summaries.append(dict(condition=condition, arm=arm, workflows=n,
                    false_splits=errors if null else None,
                    false_split_upper_bound=upper if null else None,
                    exact_null_excess_p=float(probability) if null else None,
                    candidate_rejected=bool(null and arm == 'candidate' and probability <= Fraction(1, 80)),
                    exact_chr12_recovery=recovery if not null else None,
                    mean_all_cell_assignment_accuracy=float(np.mean([r['assignment_correct']/r['assignment_n'] for r in selected])),
                    mean_tumour_calls=float(np.mean([r['tumour_calls'] for r in selected]))))
        write_json(summaries, args.out/'performance.json')
        write_json(dict(complete_workflows=60, condition_datasets=30, independent_base_seeds=10,
            candidate_rejected=any(r['candidate_rejected'] for r in summaries),
            target_error_control_established=False, biological_specimens=0,
            scientific_criteria_changed=False), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU paired workflow error and power with exact binomial pilot rule')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
