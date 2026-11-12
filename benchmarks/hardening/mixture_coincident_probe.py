"""A bounded algebraic counterexample for null-relative mixture normalization."""
import argparse
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.special import expit
from haplodose.clones import _event_mixture_objective
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    record = start_run(args.out, dict(script=Path(__file__), design=args.task/'DESIGN.md'), {},
        dict(scope='Known mathematical identity only, no model/power or biological claim'), mode='coincident_mixture_probe')
    try:
        O = E = np.full((17, 2), 2.)
        r = np.full((17, 2), 2.)
        rows = []
        for logit in (-6., -2., -.7, 0., .3, 2., 6.):
            pi = expit(logit)
            # Identical components give exactly the null distribution for all
            # observations, since these mathematical weights sum to one.
            assert Fraction(pi)+(1-Fraction(pi)) == 1
            value, gradient, q = _event_mixture_objective(np.array([0., 0., logit]), O, E, r, 0.)
            prototype_total = np.log1p(pi*np.expm1(np.zeros(len(O))))
            prototype_q = expit(logit+np.zeros(len(O)))
            assert np.all(prototype_total == 0) and np.all(prototype_q == pi)
            rows.append(dict(logit=logit, expected_null_relative_objective=0., actual_objective=float(value),
                expected_mixture_gradient=0., actual_mixture_gradient=float(gradient[2]),
                max_membership_prior_discrepancy=float(np.max(abs(q-pi))),
                prototype_objective=float(-prototype_total.sum()), prototype_mixture_gradient=float(-(prototype_q-pi).sum())))
        write_json(dict(records=rows, exact_identity_violations=sum(row['actual_objective'] != 0 or
            row['actual_mixture_gradient'] != 0 or row['max_membership_prior_discrepancy'] != 0 for row in rows),
            effect_scope='Floating-point null identity; scientific gain-budget veto remains in force',
            historical_log_variation_cause_isolated=False, production_changed=False), args.out/'report.json')
        finish_run(args.out, record, backend='CPU algebraic identical-component fixture')
    except BaseException as error:
        finish_run(args.out, record, error=error)
        raise


if __name__ == '__main__':
    main()
