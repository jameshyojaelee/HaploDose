"""Exact-binomial outward bounds for the frozen two-look clone pilot."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


ALPHA=Fraction(1,80)
TARGET=Fraction(1,20)


def cdf(k,n,p):
    p=Fraction(p)
    return sum(Fraction(math.comb(n,j))*p**j*(1-p)**(n-j)for j in range(k+1))


def upper(k,n,alpha):
    """Exact rational CDF bisection to adjacent binary64 outward endpoints."""
    alpha=Fraction(alpha)
    if k==n:return dict(value=1.,bracket=[1.,1.],exact=True)
    lo,hi=0.,1.
    while True:
        mid=(lo+hi)/2
        if mid==lo or mid==hi:break
        if cdf(k,n,Fraction(mid))>alpha:lo=mid
        else:hi=mid
    assert cdf(k,n,Fraction(lo))>=alpha>=cdf(k,n,Fraction(hi))
    return dict(value=hi,bracket=[lo,hi],exact=True)


def lower(k,n,alpha):
    if k==0:return dict(value=0.,bracket=[0.,0.],exact=True)
    bound=upper(n-k,n,alpha)
    # Subtraction itself must round outward; use exact rational conversion
    # followed by a one-ulp adjustment only if the rounded float is above it.
    exact=1-Fraction(bound['value'])
    value=float(exact)
    if Fraction(value)>exact:value=math.nextafter(value,-math.inf)
    return dict(value=value,complement_upper=bound,exact=True)


def checks():
    u=upper(0,1,ALPHA)
    assert Fraction(u['bracket'][0])<=Fraction(79,80)<=Fraction(u['bracket'][1])
    l=lower(1,1,Fraction(1,40))
    assert Fraction(l['value'])<=Fraction(1,40)
    # Replication claims established by exact integer probabilities, not a
    # rounded beta inverse or a convenient numerical tolerance.
    assert cdf(0,58,TARGET)>Fraction(1,20)>=cdf(0,59,TARGET)
    minimum=1
    while cdf(0,minimum,TARGET)>ALPHA:minimum+=1
    assert minimum==86
    assert 1-cdf(2,10,TARGET)<=ALPHA<1-cdf(1,10,TARGET)
    assert cdf(0,100,TARGET)<=ALPHA<cdf(1,100,TARGET)
    return dict(exact_rational_confidence_inversion=True,
        single_population_zero_error_minimum95=59,two_null_two_look_zero_error_minimum=minimum,
        pilot_excess_boundary_errors=3,one_error_final100_futile=True,
        n1_upper=u,n1_lower=l,biological_error_control_established=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--numerics',action='store_true')
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    inputs=dict(script=Path(__file__),design=args.task/'PILOT_DESIGN.md')
    if not args.numerics:
        for condition in('matched','state','two_clone'):
            for start in(0,5):
                inputs[f'{condition}{start}']=args.task/f'pilot_{condition}_{start:02d}_v1/results.json'
    rec=start_run(args.out,inputs,dict(error_family='2 nulls x2 looks',one_sided_alpha=str(ALPHA),
        target=str(TARGET),planned_final_n=100),
        dict(scope='Conditional simulation error/power; biological n=0'),mode='pooled_prevalence_pilot_summary')
    try:
        if args.numerics:
            write_json(checks(),args.out/'report.json')
        else:
            results=[]
            for key,path in inputs.items():
                if key not in('script','design'):results.extend(json.loads(path.read_text()))
            assert len(results)==60
            summaries,paired=[],[]
            for condition in('matched','state','two_clone'):
                by_arm={}
                for arm in('production','pooled_prevalence'):
                    rows=sorted((r for r in results if r['condition']==condition and r['arm']==arm),
                                key=lambda r:r['replicate'])
                    assert [r['replicate']for r in rows]==list(range(10))
                    assert all(r['phase']==8 and r['assignment_n']==600 for r in rows)
                    by_arm[arm]=rows
                    n=len(rows);errors=sum(r['null_error']for r in rows)
                    recovery=sum(r['exact_chr12_recovery']for r in rows)
                    null=condition!='two_clone'
                    excess=1-cdf(errors-1,n,TARGET)if errors else Fraction(1)
                    bound=upper(errors,n,ALPHA)
                    final=upper(errors,100,ALPHA)
                    possible=cdf(errors,100,TARGET)<=ALPHA
                    accuracies=[r['assignment_correct']/600 for r in rows]
                    summaries.append(dict(condition=condition,arm=arm,workflows=n,
                        false_splits=errors if null else None,
                        simultaneous_false_split_upper=bound['value']if null else None,
                        simultaneous_false_split_upper_bracket=bound['bracket']if null else None,
                        exact_excess_p=str(excess)if null else None,
                        excess_rejected=bool(null and excess<=ALPHA),
                        best_final100_upper=final['value']if null else None,
                        final100_possible=bool(possible)if null else None,
                        exact_chr12_recovery=recovery if not null else None,
                        descriptive_power_exact95=[lower(recovery,n,Fraction(1,40))['value'],
                            upper(recovery,n,Fraction(1,40))['value']]if not null else None,
                        all_cell_assignment_accuracy=accuracies,
                        mean_all_cell_assignment_accuracy=math.fsum(accuracies)/n,
                        tumour_calls=[r['tumour_calls']for r in rows],
                        zero_reported_groups=sum(r['reported_groups']==0 for r in rows)))
                if condition=='two_clone':
                    pairs=list(zip(by_arm['pooled_prevalence'],by_arm['production']))
                    assert all(a['seed']==b['seed']for a,b in pairs)
                    paired.append(dict(comparator='production',
                        candidate_only=sum(a['exact_chr12_recovery']and not b['exact_chr12_recovery']for a,b in pairs),
                        comparator_only=sum(b['exact_chr12_recovery']and not a['exact_chr12_recovery']for a,b in pairs),
                        both=sum(a['exact_chr12_recovery']and b['exact_chr12_recovery']for a,b in pairs)))
            primary=[s for s in summaries if s['arm']=='pooled_prevalence'and s['condition']!='two_clone']
            write_json(summaries,args.out/'performance.json')
            write_json(dict(complete_workflows=60,generated_datasets=30,independent_base_seeds=10,
                candidate_rejected=any(s['excess_rejected']for s in primary),
                final100_certification_possible=all(s['final100_possible']for s in primary),
                target_error_control_established=False,paired_power=paired,biological_specimens=0,
                candidate_frozen=True,scientific_criteria_changed=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU exact-rational binomial bounds and paired whole-workflow power')
    except BaseException as error:
        finish_run(args.out,rec,error=error);raise


if __name__=='__main__':main()
