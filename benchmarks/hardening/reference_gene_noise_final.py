"""Final frozen 100-replicate paired workflow error/power experiment."""
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
    inputs = dict(script=Path(__file__), design=args.task/'PILOT_DESIGN.md',
                  continuation=args.task/'FINAL_DESIGN.md')
    for condition in ('matched', 'state', 'two_clone'):
        for start in (0, 5, 10, 25, 40, 55, 70, 85):
            prefix = 'pilot' if start < 10 else 'final'
            inputs[f'{condition}{start}'] = args.task/f'{prefix}_{condition}_{start:02d}_v1/results.json'
    rec = start_run(args.out, inputs, dict(error_family='2 nulls x2 looks', alpha_per_look=.0125),
        dict(scope='Conditional simulation experiment; biological specimens=0'), mode='reference_gene_noise_final')
    try:
        results = []
        for name, path in inputs.items():
            if name not in ('script', 'design', 'continuation'):
                results.extend(json.loads(path.read_text()))
        assert len(results) == 600
        summaries = []
        for condition in ('matched', 'state', 'two_clone'):
            for arm in ('production', 'candidate'):
                selected = [r for r in results if r['condition'] == condition and r['arm'] == arm]
                assert len(selected) == 100 and sorted(r['replicate'] for r in selected) == list(range(100))
                n = len(selected)
                errors = sum(r['null_error'] for r in selected)
                null = condition != 'two_clone'
                probability = Fraction(sum(math.comb(n, j)*19**(n-j) for j in range(errors, n+1)), 20**n)
                upper = 1. if errors == n else float(beta.ppf(1-.0125, errors+1, n-errors))
                recovery = sum(r['exact_chr12_recovery'] for r in selected)
                summaries.append(dict(condition=condition, arm=arm, workflows=n,
                    false_splits=errors if null else None,
                    false_split_upper_bound=upper if null else None,
                    exact_null_excess_p=float(probability) if null else None,
                    candidate_rejected=bool(null and arm == 'candidate' and probability <= Fraction(1, 80)),
                    target_bound_passed=bool(upper <= .05) if null else None,
                    exact_chr12_recovery=recovery if not null else None,
                    recovery_lower_95=float(beta.ppf(.025, recovery, n-recovery+1)) if recovery and not null else None,
                    recovery_upper_95=(1. if recovery == n else float(beta.ppf(.975, recovery+1, n-recovery))) if not null else None,
                    mean_all_cell_assignment_accuracy=float(np.mean([r['assignment_correct']/r['assignment_n'] for r in selected])),
                    mean_tumour_calls=float(np.mean([r['tumour_calls'] for r in selected]))))
        paired = {}
        for condition in ('matched', 'state', 'two_clone'):
            by_arm = {arm: {r['replicate']: r for r in results if r['condition'] == condition and r['arm'] == arm}
                      for arm in ('production', 'candidate')}
            paired[condition] = dict(
                candidate_only_recovery=sum(by_arm['candidate'][i]['exact_chr12_recovery'] and not by_arm['production'][i]['exact_chr12_recovery'] for i in range(100)),
                production_only_recovery=sum(by_arm['production'][i]['exact_chr12_recovery'] and not by_arm['candidate'][i]['exact_chr12_recovery'] for i in range(100)))
        write_json(summaries, args.out/'performance.json')
        write_json(dict(complete_workflows=600, condition_datasets=300, independent_base_seeds=100, paired=paired,
            candidate_rejected=any(r['candidate_rejected'] for r in summaries),
            conditional_null_bounds_passed=all(r['target_bound_passed'] for r in summaries if r['condition'] != 'two_clone' and r['arm'] == 'candidate'),
            biological_error_control_established=False, biological_specimens=0,
            scientific_criteria_changed=False), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU exact binomial final error and simultaneous power report')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
