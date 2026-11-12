"""Paired patient performance effects, retaining the complete comparator family."""
import argparse
from fractions import Fraction
import math
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    paths = {arm: args.task/f'{folder}/patient_scores.tsv' for arm, folder in
        [('single', 'profile_v1'), ('ensemble', 'ensemble_v1')]}
    rec = start_run(args.out, dict(script=Path(__file__), design=args.task/'EFFECTS_DESIGN.md', **paths),
        dict(family=96, unit='patient', estimand='median paired endpoint difference'),
        dict(scope='Post-hoc exposed-development performance characterization, not independent superiority'), mode='performance_effects')
    try:
        rows, probabilities = [], []
        for arm, path in paths.items():
            scores = pd.read_csv(path, sep='\t')
            scores = scores[scores.split.eq('evaluation')]
            candidate = scores[scores.method.eq('continuous_mdl')].set_index('sample')
            assert candidate.index.is_unique and len(candidate) == 24
            controls = sorted(set(scores.method)-{'continuous_mdl'})
            assert len(controls) == 12
            for control in controls:
                comparator = scores[scores.method.eq(control)].set_index('sample').loc[candidate.index]
                assert comparator.index.is_unique
                for endpoint in ('pearson', 'recall_loss', 'fpr_loss', 'precision_loss'):
                    a, b = candidate[endpoint].to_numpy(), comparator[endpoint].to_numpy()
                    ok = np.isfinite(a) & np.isfinite(b)
                    differences = (a-b if endpoint != 'fpr_loss' else b-a)[ok]
                    n, positive = len(differences), int(sum(differences > 0))
                    probability = Fraction(sum(math.comb(n, j) for j in range(positive, n+1)), 2**n)
                    probabilities.append(probability)
                    ordered = np.sort(differences)
                    indices = [j for j in range(n) if Fraction(sum(math.comb(n, k) for k in range(j+1)), 2**n) <= Fraction(1, 20*96)]
                    lower = float(ordered[max(indices)]) if indices else None
                    rows.append(dict(arm=arm, comparator=control, endpoint=endpoint,
                        paired_patients=n, excluded_patients=candidate.index[~ok].tolist(),
                        mean_improvement=float(np.mean(differences)) if n else None,
                        median_improvement=float(np.median(differences)) if n else None,
                        positive_patients=positive, exact_one_sided_sign_p=float(probability),
                        simultaneous_median_lower_bound=lower,
                        bound_status='unbounded_below' if lower is None else 'finite',
                        interpretation='conditional on paired defined endpoints; precision exclusions explicit'))
        assert len(rows) == 96
        running = Fraction(0)
        for rank, index in enumerate(sorted(range(96), key=lambda i: probabilities[i])):
            running = max(running, min(Fraction(1), probabilities[index]*(96-rank)))
            rows[index]['holm_p'] = float(running)
            rows[index]['median_superiority_descriptive'] = running <= Fraction(1, 20)
        write_json(rows, args.out/'effects.json')
        pd.DataFrame(rows).to_csv(args.out/'effects.tsv', sep='\t', index=False)
        write_json(dict(contrasts=96, patient_unit=True, comparator_family_complete=True,
            descriptive_median_superiority=sum(r['median_superiority_descriptive'] for r in rows),
            prospective_claim=False, scientific_criteria_changed=False), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU paired patients, exact integer binomial tails and Holm family')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
