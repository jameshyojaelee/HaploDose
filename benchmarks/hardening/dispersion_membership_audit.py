"""Exact rational NB-mixture moments, independent of SciPy probability evaluation."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.provenance import start_run, finish_run


def exact(means, prior, N, subtraction='observed'):
    r = 5
    probabilities = [(Fraction(r, r+mu))**r for mu in means]
    softnum, softden, hardnum, hardden = [Fraction(0) for _ in range(4)]
    assignments = []
    for y in range(N+1):
        joint = [pi*p for pi, p in zip(prior, probabilities)]
        moment = [(y-mu)**2-(y if subtraction == 'observed' else mu) for mu in means]
        # Mixture probability times posterior responsibility cancels exactly.
        softnum += sum(p*m for p, m in zip(joint, moment))
        softden += sum(p*mu**2 for p, mu in zip(joint, means))
        selected = int(joint[1] > joint[0])
        assignments.append(selected)
        hardnum += sum(joint)*moment[selected]
        hardden += sum(joint)*means[selected]**2
        probabilities = [p*Fraction(y+r, y+1)*Fraction(mu, r+mu) for p, mu in zip(probabilities, means)]
    return dict(soft=(softnum, softden), MAP=(hardnum, hardden)), np.array(assignments)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--panel-name', default='panel_v1')
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    panel = a.task/a.panel_name
    run = start_run(a.out, dict(script=Path(__file__), design=a.task/'DESIGN.md',
        moments=panel/'moments.tsv', report=panel/'report.json', assignment=panel/'assignments.npz'), {},
        dict(scope='Exact rational moment check; no biological or workflow validation'), mode='dispersion_membership_audit')
    try:
        report = json.loads((panel/'report.json').read_text())
        frame = pd.read_csv(panel/'moments.tsv', sep='\t')
        saved = np.load(panel/'assignments.npz')
        assert len(frame) == 18 and len(report['cases']) == 9
        rows = []
        for case in report['cases']:
            means = [int(v) for v in case['means']]
            numerator_prior = int(case['case'].split('_p')[1])
            prior = [Fraction(numerator_prior, 10), Fraction(10-numerator_prior, 10)]
            result, assignments = exact(means, prior, case['N'], report.get('subtraction', 'observed'))
            assert np.array_equal(assignments, saved[case['case']+'__assignment'])
            for method, (num, den) in result.items():
                row = frame[(frame.case == case['case']) & (frame.method == method)].iloc[0]
                estimate = float(num/den)
                discrepancy = abs(estimate-row.estimated_phi)
                assert discrepancy <= row.numerical_budget
                assert abs(float(num)-row.numerator) <= row.numerical_budget*max(1., abs(float(num)))
                assert abs(float(den)-row.denominator) <= row.numerical_budget*max(1., abs(float(den)))
                if method == 'soft':
                    assert abs(estimate-.2) <= row.ratio_tail_bound+row.numerical_budget
                rows.append(dict(case=case['case'], method=method, exact_truncated_phi=estimate,
                    discrepancy=discrepancy, true_phi=.2, signed_error=estimate-.2,
                    numerical_budget=float(row.numerical_budget), tail_bound=float(row.ratio_tail_bound)))
        pd.DataFrame(rows).to_csv(a.out/'exact_moments.tsv', sep='\t', index=False)
        (a.out/'report.json').write_text(json.dumps(dict(cases=9, exact_ratios=18,
            exact_MAP_decisions=sum(c['N']+1 for c in report['cases']), biological_specimens=0,
            all_checks_passed=True, production_changed=False, scientific_criteria_changed=False), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, run, error=error)
        raise
    finish_run(a.out, run, backend='CPU Python Fraction NB recurrence and exact membership decisions')


if __name__ == '__main__':
    main()
