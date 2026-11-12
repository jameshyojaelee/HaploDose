"""Prespecified stricter rank operating points; retain every prior endpoint.

This compares additional candidate abstention rules on a development simulation
panel. It never changes the original alpha=.05 run or selects a production rule.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq,linear_sum_assignment
from scipy.stats import beta,binom
from haplodose.provenance import start_run,finish_run,file_record

FAMILY=72  # original54 and3strictercutoffs x2banks x3conditions
ALPHAS=(.005,.01,.025)


def confidence(k,n,power):
    tail=.05/FAMILY
    if power:
        value=0. if k==0 else float(beta.ppf(tail,k,n-k+1))
        reference=0. if k==0 else brentq(lambda p:binom.sf(k-1,n,p)-tail,0.,1.,xtol=1e-14)
    else:
        value=1. if k==n else float(beta.ppf(1-tail,k+1,n-k))
        reference=1. if k==n else brentq(lambda p:binom.cdf(k,n,p)-tail,0.,1.,xtol=1e-14)
    # Root solver absolute tolerance plus inverse-CDF rounding; recorded as a
    # numerical reference discrepancy, never a scientific acceptance tolerance.
    if abs(value-reference)>1e-14+16*np.finfo(float).eps:raise ArithmeticError('Independent binomial inverse differs')
    return value,reference


def membership(record):
    artifact=record['assignments']
    if file_record(artifact['path'])['sha256']!=artifact['sha256']:raise ValueError('Assignment artifact changed')
    with np.load(artifact['path']) as z:
        truth=z['simulated_CNA_labels'];labels=z['raw_MAP'];eligible=z['eligible_tumour_groups']
    if len(truth)!=600:raise ValueError('Fixed specimen size changed')
    counts=np.array([[np.count_nonzero((truth==a)&(labels==b)) for b in eligible] for a in np.unique(truth)])
    if counts.shape[1]:
        a,b=linear_sum_assignment(-counts);accuracy=float(counts[a,b].sum()/len(truth))
    else:accuracy=0.
    return dict(accuracy=accuracy,coverage=float(np.isin(labels,eligible).mean()),groups=len(eligible))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--source-task',type=Path,required=True);args=parser.parse_args()
    source=args.source_task
    audit=json.loads((source/'audit_v1/report.json').read_text())
    inputs={'script':__file__,'audit':source/'audit_v1/report.json',
            'summary':source/'summary_v1/summary.json','evaluation':source/'summary_v1/evaluation.json'}
    rows={}
    for role,n in [('cal_nb',199),('cal_state',199),('eval_nb',200),('eval_state',200),('eval_two',200)]:
        directory=source/'shards_v1'/role
        inputs[role]=directory/'records.jsonl';inputs[role+'_manifest']=directory/'run_manifest.json'
        if json.loads(inputs[role+'_manifest'].read_text())['status']!='complete':raise ValueError('Incomplete shard')
        rows[role]=[json.loads(line) for line in inputs[role].read_text().splitlines()]
        if len(rows[role])!=n:raise ValueError('Specimen population changed')
    args.out.mkdir(exist_ok=False)
    rec=start_run(args.out,inputs,dict(family=FAMILY,additional_alphas=ALPHAS,original_alpha=.05),
        dict(scope='Development stricter operating-point comparison; no rule selected or promoted',
             scientific_target='Whole-specimen false-split upper bound<=.05, with simultaneous detection and membership assessment'),
        mode='pooled_threshold_tradeoff')
    try:
        old=json.loads(inputs['summary'].read_text());retained=[]
        for row in old['summary']+old['retained_historical']:
            value,reference=confidence(row['events'],row['n'],row['direction']=='lower')
            retained.append(dict(row,simultaneous_bound=value,independent_bound=reference,family=FAMILY))
        if len(retained)!=54:raise ValueError('Prior confidence family incomplete')
        old_decisions=json.loads(inputs['evaluation'].read_text());summary=[];decisions=[]
        memberships={r['replicate']:membership(r) for r in rows['eval_two']}
        for role,condition in [('eval_nb','matched_nb'),('eval_state','correlated_state'),('eval_two','two_clones')]:
            for bank_name,method in [('cal_nb','NB_calibrated_gate'),('cal_state','state_calibrated_gate')]:
                bank=np.array([r['score'] for r in rows[bank_name]])
                old_rows={r['replicate']:r for r in old_decisions if r['condition']==condition and r['method']==method}
                if len(old_rows)!=200:raise ValueError('Original decisions incomplete')
                for alpha in ALPHAS:
                    accepted=[]
                    for row in rows[role]:
                        p=float((1+np.count_nonzero(bank>=row['score']))/200)
                        valid=row['raw_candidate_count']>1 and row['score']>0
                        if p!=old_rows[row['replicate']]['rank_p'] or bool(valid and p<=.05)!=old_rows[row['replicate']]['event']:
                            raise ValueError('Original operating decisions changed')
                        event=bool(valid and p<=alpha);accepted.append(event)
                        decisions.append(dict(condition=condition,method=method,alpha=alpha,
                            replicate=row['replicate'],event=event,rank_p=p,raw_groups=row['raw_candidate_count']))
                    power=condition=='two_clones';k=sum(accepted);value,reference=confidence(k,200,power)
                    item=dict(condition=condition,method=method,alpha=alpha,events=k,n=200,
                        direction='lower' if power else 'upper',simultaneous_bound=value,independent_bound=reference,
                        family=FAMILY,error_criterion_met=None if power else value<=.05,
                        original_alpha_events=sum(r['event'] for r in old_rows.values()))
                    if power:
                        item.update(exactly_two_and_accepted=sum(flag and memberships[r['replicate']]['groups']==2
                                                              for flag,r in zip(accepted,rows[role])),
                            mean_coverage_if_gate_applied=float(np.mean([memberships[r['replicate']]['coverage'] if flag else 0.
                                                              for flag,r in zip(accepted,rows[role])])),
                            mean_correct_fraction_all_cells_if_gate_applied=float(np.mean([memberships[r['replicate']]['accuracy'] if flag else 0.
                                                              for flag,r in zip(accepted,rows[role])])))
                    summary.append(item)
        (args.out/'decisions.json').write_text(json.dumps(decisions,indent=2))
        (args.out/'summary.json').write_text(json.dumps(dict(summary=summary,retained_historical=retained,
            family=FAMILY,original_criterion=.05,selection='No default selected; stricter operating rules do not repair memberships',
            limits='Known-reference expression-only simulation; law matching required; no gene/allele/spatial/raw-preprocessing or biological guarantee'),indent=2))
        finish_run(args.out,rec,backend='cpu')
    except Exception as error:
        finish_run(args.out,rec,error=error);raise


if __name__=='__main__':main()
