"""Monte Carlo precision of already-frozen bootstrap reports, not new patient inference."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from scipy.stats import beta, binom
from haplodose.provenance import start_run, finish_run, file_record


def holm(values):
    return [min(1., max((len(values)-sum(v < u for v in values))*u
        for u in values if u <= p)) for p in values]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--previous',type=Path,required=True)
    p.add_argument('--audit',type=Path,required=True)
    p.add_argument('--design',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    audit_job=a.audit.parent/'analysis_22972645'
    inputs=dict(script=Path(__file__),summary=a.previous/'summary.json',
        previous_manifest=a.previous/'run_manifest.json',audit=a.audit/'report.json',
        audit_job=audit_job/'job.tsv',audit_execution=audit_job/'execution.json',
        audit_stdout=audit_job/'output.txt',design=a.design)
    m=json.loads(inputs['previous_manifest'].read_text())
    if m['status']!='complete' or file_record(inputs['summary'])['sha256']!=m['output']['files']['summary.json']['sha256']:
        raise AssertionError('Prior completed evidence changed')
    job=dict(line.split('\t',1) for line in inputs['audit_job'].read_text().splitlines())
    execution=json.loads(inputs['audit_execution'].read_text())
    if job['exit']!='0' or Path(job['output']).resolve()!=a.audit.resolve():
        raise AssertionError('Historical audit did not complete for the requested report')
    if job['script']!='benchmarks/hardening/gene_aggregation_audit.py' or execution['source']!=job['source']:
        raise AssertionError('Historical audit execution identity differs')
    if json.loads(inputs['audit'].read_text())!=json.loads(inputs['audit_stdout'].read_text()):
        raise AssertionError('Historical audit report differs from captured output')
    rec=start_run(a.out,inputs,dict(scientific_family=93,mc_family=279,mc_confidence=.95),
        dict(scope='Conditional Monte Carlo precision only; no new biological evidence or modified inference'),
        mode='profile_mc_precision')
    try:
        old=json.loads(inputs['summary'].read_text());rows=old['comparisons']
        if len(rows)!=93 or old['params']['bootstrap']!=20000:
            raise AssertionError('Unexpected frozen bootstrap population')
        # 93 bootstrap exceedance probabilities and two quantiles per interval.
        tail=.05/(2*279);result=[];nominal=[];lower=[];upper=[]
        for row in rows:
            n=int(row['valid_bootstrap_replicates'])
            if not row['testable']:
                nominal.append(1.);lower.append(1.);upper.append(1.)
                result.append(dict(endpoint=row['endpoint'],a=row['a'],b=row['b'],testable=False,
                    reason='Original untestable result retained'))
                continue
            k=round(row['p']*(n+1)-1)
            if not 0<=k<=n or float(Fraction(k+1,n+1))!=row['p']:
                raise AssertionError('Cannot reconstruct exact recorded exceedance count')
            lo=0. if k==0 else float(beta.ppf(tail,k,n-k+1))
            hi=1. if k==n else float(beta.isf(tail,k+1,n-k))
            quantiles=[]
            for q in (.05/(2*93),1-.05/(2*93)):
                # X_(0)=-infinity and X_(n+1)=+infinity. Binomial order-statistic
                # coverage is exact for continuous laws, conservative with ties.
                left=int(binom.ppf(tail,n,q));right=int(binom.isf(tail,n,q))+1
                if (left>0 and binom.cdf(left-1,n,q)>tail) or (right<=n and binom.sf(right-1,n,q)>tail):
                    raise AssertionError('Order-statistic coverage rank is invalid')
                quantiles.append(dict(probability=q,expected_tail_draws=n*min(q,1-q),
                    lower_one_based_rank=left,upper_one_based_rank=right,
                    finite_lower=left>0,finite_upper=right<=n))
            result.append(dict(endpoint=row['endpoint'],a=row['a'],b=row['b'],testable=True,n=n,
                exceedances=k,recorded_p=row['p'],recorded_holm=row['p_holm'],
                bootstrap_exceedance_probability_mc_lower=lo,bootstrap_exceedance_probability_mc_upper=hi,
                interval_quantile_mc_ranks=quantiles,
                zero_count_upper_analytic=None if k else -math.expm1(math.log(tail)/n),
                recorded_superior=row['superior']))
            nominal.append(row['p']);lower.append(lo);upper.append(hi)
        adjusted_lower,adjusted_upper=holm(lower),holm(upper)
        for row,lo,hi in zip(result,adjusted_lower,adjusted_upper):
            row.update(holm_mc_lower=lo,holm_mc_upper=hi)
        (a.out/'report.json').write_text(json.dumps(dict(results=result,mc_family=279,
            scientific_family=93,mc_confidence=.95,biological_units_added=0,
            interpretation='Numerical uncertainty conditional on the original data and eight-reference ensemble; '
                'does not validate bootstrap coverage, reference adequacy, DNA truth or a biological confidence claim',
            original_results_changed=False),indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
