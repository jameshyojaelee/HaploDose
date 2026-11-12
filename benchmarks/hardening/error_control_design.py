"""Outcome-free precision audit of the frozen specimen error criterion.

No observed simulation or biological results are read. This does not amend the
running experiment, choose a new cutoff, or create independent replication.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.stats import beta,binom
from haplodose.provenance import start_run,finish_run


def upper(k,n,alpha):return 1. if k==n else float(beta.ppf(1-alpha,k+1,n-k))


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    a.out.mkdir(exist_ok=False)
    parameters=dict(target=.05,evaluation_n=200,calibration_n=199,families=[1,27,45])
    rec=start_run(a.out,{'source':__file__},parameters,
        dict(scope='Exact binomial/order-statistic arithmetic only; no outcomes, changed operating point, or biological evidence'),mode='error_control_precision_design')
    try:
        rows=[];target=.05;n=200;m=199
        for family in parameters['families']:
            alpha=.05/family
            acceptable=[k for k in range(n+1) if upper(k,n,alpha)<=target]
            critical=max(acceptable) if acceptable else -1
            if critical<0:raise ArithmeticError('No acceptable error count')
            # Independent equivalence: CP upper <= target iff the lower tail
            # at target is <= alpha, with the discrete boundary retained.
            tail_critical=max(k for k in range(n+1) if binom.cdf(k,n,target)<=alpha)
            if critical!=tail_critical:raise ArithmeticError('Beta/binomial inversion disagreement')
            reference=brentq(lambda x:binom.cdf(critical,n,x)-alpha,0.,1.,
                xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps)
            discrepancy=abs(reference-upper(critical,n,alpha))
            if discrepancy>64*np.finfo(float).eps:raise ArithmeticError('Independent numerical bound failed')
            # These are prospective operating-rate scenarios, never estimates
            # from the ongoing fits. Include the boundary and fractions of it.
            scenarios=[dict(true_error_rate=float(rate),probability_of_passing=float(binom.cdf(critical,n,rate)))
                       for rate in (target,target/2,target/5,target/10)]
            design_rates={}
            for desired in (.8,.9):
                rate=brentq(lambda x:binom.cdf(critical,n,x)-desired,0.,target,
                            xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps)
                design_rates[str(desired)]=rate
            # For continuous exchangeable null scores, the k-th largest of m
            # calibration scores has population exceedance mass Beta(k,m+1-k).
            # Ties + strict exceedance are conservative. This is an alternative
            # design calculation, not an amendment to the frozen rank-p rule.
            orders=[k for k in range(1,m+1) if beta.ppf(1-alpha,k,m+1-k)<=target]
            largest=max(orders) if orders else None
            rows.append(dict(family=family,tail_alpha=alpha,evaluation_n=n,maximum_acceptable_errors=critical,
                upper_at_boundary=upper(critical,n,alpha),upper_with_one_more_error=upper(critical+1,n,alpha),
                independent_bound_discrepancy=discrepancy,scenarios=scenarios,
                largest_true_rate_for_validation_pass_probability=design_rates,
                calibration_order_statistic=dict(n=m,largest_acceptable_descending_order=largest,
                    population_tail_upper=None if largest is None else float(beta.ppf(1-alpha,largest,m+1-largest)),
                    expected_tail_mass=None if largest is None else largest/(m+1),
                    assumptions='Independent exchangeable continuous null scores; ties with strict exceedance conservative; no model-misspecification protection')))
        single_min=int(np.ceil(np.log(.05)/np.log(1-target)))
        if single_min!=59 or upper(0,59,.05)>target or upper(0,58,.05)<=target:
            raise ArithmeticError('Biological minimum replication check failed')
        report=dict(parameters=parameters,calculations=rows,minimum_independent_null_specimens_zero_errors_single_endpoint=single_min,
            distinction='Failure to obtain an upper bound <=5% is not itself evidence that true error exceeds5%; a lower bound above5% would support that conclusion.',
            decision='Keep current fixed counts, alpha and scientific pass criterion unchanged. No operating-point selection from this report.',
            limits=['Independent biological specimens cannot be replaced by these calculations or simulated seeds.',
                    'Calibration-bank risk bounds need exchangeability and independent calibration; RNA-state mismatch can violate applicability.',
                    'A lower operating error can cost detection power; neither this arithmetic nor abstention establishes clone recovery.'])
        (a.out/'summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
