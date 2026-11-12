"""Prespecified two-look exact error rule and simultaneous recovery reporting."""
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
            inputs[f'{condition}{start}'] = args.task/f'pilot_{condition}_{start:02d}_v3/results.json'
    rec = start_run(args.out, inputs, dict(error_family='2 nulls x2 looks',
        one_sided_alpha=.0125, target=.05, planned_final_n=100),
        dict(scope='Conditional simulation test; not independent biological validation'),
        mode='conditioned_gene_pilot_summary')
    try:
        results = []
        for key, path in inputs.items():
            if key not in ('script', 'design'):
                results.extend(json.loads(path.read_text()))
        assert len(results) == 90
        summaries, paired = [], []
        for condition in ('matched', 'state', 'two_clone'):
            by_arm = {}
            for arm in ('production', 'depth_conditioned', 'depth_gene_score'):
                rows = sorted((r for r in results if r['condition'] == condition and r['arm'] == arm),
                    key=lambda r: r['replicate'])
                assert [r['replicate'] for r in rows] == list(range(10))
                assert all(r['phase'] == 6 and r['assignment_n'] == 600 for r in rows)
                by_arm[arm] = rows
                n = len(rows)
                errors = sum(r['null_error'] for r in rows)
                recovery = sum(r['exact_chr12_recovery'] for r in rows)
                null = condition != 'two_clone'
                tail = Fraction(sum(math.comb(n, j)*19**(n-j) for j in range(errors, n+1)), 20**n)
                upper = 1. if errors == n else float(beta.ppf(1-.0125, errors+1, n-errors))
                final_best_upper = float(beta.ppf(1-.0125, errors+1, 100-errors))
                power_ci = [0. if recovery == 0 else float(beta.ppf(.025, recovery, n-recovery+1)),
                    1. if recovery == n else float(beta.ppf(.975, recovery+1, n-recovery))]
                summaries.append(dict(condition=condition, arm=arm, workflows=n,
                    false_splits=errors if null else None,
                    simultaneous_false_split_upper=upper if null else None,
                    exact_excess_p=float(tail) if null else None,
                    excess_rejected=bool(null and tail <= Fraction(1, 80)),
                    best_final100_upper=final_best_upper if null else None,
                    final100_possible=bool(final_best_upper <= .05) if null else None,
                    exact_chr12_recovery=recovery if not null else None,
                    descriptive_power_exact95=power_ci if not null else None,
                    all_cell_assignment_accuracy=[r['assignment_correct']/600 for r in rows],
                    mean_all_cell_assignment_accuracy=float(np.mean([r['assignment_correct']/600 for r in rows])),
                    tumour_calls=[r['tumour_calls'] for r in rows],
                    zero_reported_groups=sum(r['reported_groups'] == 0 for r in rows)))
            if condition == 'two_clone':
                for comparator in ('production', 'depth_conditioned'):
                    pairs = list(zip(by_arm['depth_gene_score'], by_arm[comparator]))
                    assert all(a['seed'] == b['seed'] for a, b in pairs)
                    paired.append(dict(comparator=comparator,
                        candidate_only=sum(a['exact_chr12_recovery'] and not b['exact_chr12_recovery'] for a,b in pairs),
                        comparator_only=sum(b['exact_chr12_recovery'] and not a['exact_chr12_recovery'] for a,b in pairs),
                        both=sum(a['exact_chr12_recovery'] and b['exact_chr12_recovery'] for a,b in pairs)))
        primary = [s for s in summaries if s['arm'] == 'depth_gene_score' and s['condition'] != 'two_clone']
        write_json(summaries, args.out/'performance.json')
        write_json(dict(complete_workflows=90, generated_datasets=30, independent_base_seeds=10,
            candidate_rejected=any(s['excess_rejected'] for s in primary),
            final100_certification_possible=all(s['final100_possible'] for s in primary),
            target_error_control_established=False, paired_power=paired,
            biological_specimens=0, candidate_frozen=True, scientific_criteria_changed=False),
            args.out/'report.json')
        finish_run(args.out, rec, backend='CPU exact two-look simultaneous null bounds and paired power')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
