"""Independent recount, binomial inversion and paired regression identity audit."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.stats import binom, binomtest


def read(path):return [json.loads(line) for line in path.read_text().splitlines()]
def key(row):return row['condition'],row['role'],row['replicate']
def record(path):return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def bound(k,n,alpha,direction):
    if direction=='upper':
        if k==n:return 1.
        return brentq(lambda p:binom.cdf(k,n,p)-alpha,0.,1.,xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps)
    if k==0:return 0.
    return brentq(lambda p:binom.sf(k-1,n,p)-alpha,0.,1.,xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True);p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--predictive-dir',default='predictive_rank_v3')
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    original={key(r):r for name in ['calibration.jsonl','evaluation.jsonl'] for r in read(a.baseline/name)}
    report={'scope':'Paired inspected simulations, fixed biological claims and thresholds; no new biological replication','sources':[]}
    for name in ['rank_v1',a.predictive_dir]:
        root=a.task/name
        rows={key(r):r for name in ['calibration.jsonl','evaluation.jsonl'] for r in read(root/name)}
        if set(rows)!=set(original):raise ValueError('Specimen identities differ')
        input_mismatch=[k for k,r in rows.items() if r['input_sha256']!=original[k]['input_sha256']]
        if input_mismatch:raise ValueError('Generator/input changed')
        comparisons=[];changed_labels=0;posterior_max_difference=0.
        for k,r in rows.items():
            new=np.load(r['assignments']['path']);old=np.load(original[k]['assignments']['path'])
            for arr in ['simulated_CNA_labels','latent_state']:
                if not np.array_equal(new[arr],old[arr]):raise ValueError('Latent truth changed')
            changed_labels+=int(not np.array_equal(new['raw_MAP'],old['raw_MAP']))
            if name==a.predictive_dir:
                if new['raw_R'].shape!=old['raw_R'].shape:raise ValueError('Observer changed fitted dimensions')
                posterior_max_difference=max(posterior_max_difference,float(np.abs(new['raw_R']-old['raw_R']).max()))
        if name==a.predictive_dir and changed_labels:
            raise ValueError('Observer changed raw assignments; inspect before interpreting membership statistic')
        report[name]=dict(identical_input_specimens=len(rows),changed_raw_assignments=changed_labels,
                         maximum_probability_difference=posterior_max_difference if name==a.predictive_dir else None)
        summary=json.loads((root/'summary.json').read_text())['summary']
        evaluation=read(root/'evaluation.jsonl')
        cal=read(root/'calibration.jsonl')
        for row in evaluation:
            for method,c in [('NB_calibrated_gate','matched_nb'),('state_calibrated_gate','correlated_state')]:
                bank=[r['score'] for r in cal if r['condition']==c]
                prob=(1+sum(s>=row['score'] for s in bank))/200
                if prob!=row['methods'][method]['rank_p']:raise ValueError('Rank discrepancy')
                eligible=row['conditional_proposal_count']>0 if name==a.predictive_dir else row['raw_candidate_count']>1
                verdict=eligible and row['score']>0 and prob<=.05
                if verdict!=row['methods'][method]['accepted_multigroup']:raise ValueError('Gate discrepancy')
        for row in summary:
            rr=[r for r in evaluation if r['condition']==row['condition']]
            events=sum(r['methods'][row['method']]['accepted_multigroup'] for r in rr)
            if events!=row['events'] or len(rr)!=row['n']:raise ValueError('Count discrepancy')
            ref=bound(events,len(rr),.05/9,row['joint_family9_direction'])
            comparisons.append(dict(condition=row['condition'],method=row['method'],events=events,n=len(rr),
                independent_bound=ref,saved_bound=row['joint_family9_bound'],difference=row['joint_family9_bound']-ref))
        report[name]['family9']=comparisons
        report['sources'] += [record(root/name) for name in ['summary.json','calibration.jsonl','evaluation.jsonl']]
    root=a.task/a.predictive_dir
    comparison=json.loads((root/'membership_comparison.json').read_text());rows=read(root/'comparison_evaluation.jsonl')
    cal=read(root/'calibration.jsonl')
    for row in rows:
        for method,c in [('hard_NB_calibrated_gate','matched_nb'),('hard_state_calibrated_gate','correlated_state')]:
            prob=(1+sum(r['hard_score']>=row['hard_score'] for r in cal if r['condition']==c))/200
            if row['methods'][method]['rank_p']!=prob:raise ValueError('Hard rank discrepancy')
            expected=row['conditional_proposal_count']>0 and row['hard_score']>0 and prob<=.05
            if row['methods'][method]['accepted_multigroup']!=expected:raise ValueError('Hard gate discrepancy')
    family=[]
    for row in comparison['summary']:
        rr=[r for r in rows if r['condition']==row['condition']]
        count=sum(r['methods'][row['method']]['accepted_multigroup'] for r in rr)
        if count!=row['events'] or len(rr)!=row['n']:raise ValueError('Family15 count discrepancy')
        b=bound(count,len(rr),.05/15,row['direction'])
        family.append(dict(**row,independent_bound=b,difference=row['simultaneous_family15_bound']-b))
    paired=[]
    for condition in comparison['design']['conditions']:
        rr=[r for r in rows if r['condition']==condition]
        for method in ['NB_calibrated_gate','state_calibrated_gate']:
            soft=np.array([r['methods'][method]['accepted_multigroup'] for r in rr])
            hard=np.array([r['methods']['hard_'+method]['accepted_multigroup'] for r in rr])
            soft_only=int((soft & ~hard).sum());hard_only=int((hard & ~soft).sum());n=soft_only+hard_only
            paired.append(dict(condition=condition,method=method,soft_only=soft_only,hard_only=hard_only,
                difference=(soft_only-hard_only)/len(rr),nominal_p=1. if n==0 else float(binomtest(soft_only,n,.5).pvalue)))
    running=0.
    for rank,i in enumerate(sorted(range(len(paired)),key=lambda i:paired[i]['nominal_p'])):
        running=max(running,min(1.,(len(paired)-rank)*paired[i]['nominal_p']))
        paired[i]['holm_family6_p']=running
    report['family15']=family;report['paired_comparison']=paired
    report['sources'] += [record(root/name) for name in ['membership_comparison.json','comparison_evaluation.jsonl']]
    report['status']='completed_independent_recount'
    (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],input_specimens_per_experiment=998,paired_comparison=paired)))


if __name__=='__main__':main()
