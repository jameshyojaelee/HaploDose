"""Fixed99-endpoint summary and independent audit of the corrected-wheel replay."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq,linear_sum_assignment
from scipy.stats import beta,binom,binomtest
from haplodose.provenance import start_run,finish_run,file_record

ROLES={'cal_nb':('matched_nb',199),'cal_state':('correlated_state',199),
       'eval_nb':('matched_nb',200),'eval_state':('correlated_state',200),'eval_two':('two_clones',200)}
RULES=[('raw',None,None)]+[(method,bank,alpha) for method,bank in
       [('NB_calibrated_gate','cal_nb'),('state_calibrated_gate','cal_state')] for alpha in (.05,.005,.01,.025)]
FAMILY=99


def read(path):return [json.loads(line) for line in path.read_text().splitlines()]


def verify(record):
    if file_record(record['path'])['sha256']!=record['sha256']:raise AssertionError('Frozen artifact changed')


def collect(a):
    inputs=dict(script=__file__,old_summary=a.old_tradeoff/'tradeoff_v1/summary.json',
                old_strict_decisions=a.old_tradeoff/'tradeoff_v1/decisions.json',
                old_decisions=a.old_task/'summary_v1/evaluation.json',old_audit=a.old_task/'audit_v1/report.json')
    rows={}
    for role,(_,n) in ROLES.items():
        rows[role]=[]
        for start in (0,50,100,150):
            directory=a.task/'chunks_v1'/f'{role}_{start}'
            manifest_path=directory/'run_manifest.json';manifest=json.loads(manifest_path.read_text())
            if manifest['status']!='complete':raise ValueError('Incomplete chunk')
            params=manifest['execution']['parameters']
            if (params['role'],params['start'],params['stop'],params['confidence_family'])!=(role,start,min(start+50,n),99):
                raise ValueError('Chunk design differs')
            path=directory/'records.jsonl';part=read(path)
            if [r['replicate'] for r in part]!=list(range(start,min(start+50,n))):raise ValueError('Missing/duplicate/reordered specimen')
            rows[role].extend(part);inputs[f'{role}:{start}:records']=path;inputs[f'{role}:{start}:manifest']=manifest_path
        inputs[f'old:{role}']=a.old_task/'shards_v1'/role/'records.jsonl'
    return rows,inputs


def quantile(k,n,power):
    tail=.05/FAMILY
    return (0. if k==0 else float(beta.ppf(tail,k,n-k+1))) if power else (
        1. if k==n else float(beta.ppf(1-tail,k+1,n-k)))


def summaries(a,rows,inputs):
    old=json.loads(inputs['old_summary'].read_text())
    if old['family']!=72:raise ValueError('Previous family changed')
    retained=[]
    for r in old['retained_historical']+old['summary']:
        value=quantile(r['events'],r['n'],r['direction']=='lower')
        # Family widening changes both the bound and its decision. Do not carry
        # the earlier independent inverse-CDF value as if it checked this bound.
        updated={k:v for k,v in r.items() if k!='independent_bound'}
        updated.update(simultaneous_bound=value,family=FAMILY,
            error_criterion_met=None if r['direction']=='lower' else value<=.05,
            historical_bound=dict(family=r['family'],simultaneous_bound=r['simultaneous_bound'],
                independent_bound=r.get('independent_bound'),error_criterion_met=r.get('error_criterion_met')))
        retained.append(updated)
    if len(retained)!=72:raise ValueError('Previous endpoint omitted')
    prior=json.loads(inputs['old_strict_decisions'].read_text())
    prior += [dict(r,alpha=None if r['method']=='raw' else .05) for r in json.loads(inputs['old_decisions'].read_text())]
    old_index={(r['condition'],r['method'],r['alpha'],r['replicate']):r['event'] for r in prior}
    if len(old_index)!=5400:raise ValueError('Old paired operating points incomplete')
    membership={}
    for row in rows['eval_two']:
        verify(row['assignments'])
        with np.load(row['assignments']['path']) as z:
            truth=z['simulated_CNA_labels'];label=z['raw_MAP'];eligible=z['eligible_tumour_groups']
        count=np.array([[np.count_nonzero((truth==t)&(label==g)) for g in eligible] for t in np.unique(truth)])
        if len(eligible):
            x,y=linear_sum_assignment(-count);accuracy=float(count[x,y].sum()/len(truth))
        else:accuracy=0.
        membership[row['replicate']]=dict(accuracy=accuracy,coverage=float(np.isin(label,eligible).mean()),groups=len(eligible))
    summary=[];decisions=[];paired=[]
    for role in ('eval_nb','eval_state','eval_two'):
        condition,n=ROLES[role];power=role=='eval_two'
        for method,bank_name,alpha in RULES:
            bank=None if bank_name is None else np.asarray([r['score'] for r in rows[bank_name]])
            events=[];lost=gained=0
            for r in rows[role]:
                probability=None if bank is None else float((1+np.count_nonzero(bank>=r['score']))/200)
                event=bool(r['raw_candidate_count']>1 and (bank is None or r['score']>0 and probability<=alpha))
                before=old_index[(condition,method,alpha,r['replicate'])]
                lost+=int(before and not event);gained+=int(event and not before);events.append(event)
                decisions.append(dict(stage='precision_pooled',condition=condition,method=method,alpha=alpha,
                    replicate=r['replicate'],event=event,rank_p=probability,score=r['score'],raw_groups=r['raw_candidate_count']))
            k=sum(events);value=quantile(k,n,power)
            row=dict(stage='precision_pooled',condition=condition,method=method,alpha=alpha,events=k,n=n,
                direction='lower' if power else 'upper',simultaneous_bound=value,family=FAMILY,
                error_criterion_met=None if power else value<=.05,
                raw_exactly_two=sum(r['raw_candidate_count']==2 for r in rows[role]),
                raw_ARI_mean=float(np.mean([r['ARI_against_simulated_CNA_labels'] for r in rows[role]])))
            if power:
                row.update(exactly_two_and_accepted=sum(flag and membership[r['replicate']]['groups']==2 for flag,r in zip(events,rows[role])),
                    mean_correct_fraction_all_cells_if_gate_applied=float(np.mean([membership[r['replicate']]['accuracy'] if flag else 0. for flag,r in zip(events,rows[role])])),
                    mean_coverage_if_gate_applied=float(np.mean([membership[r['replicate']]['coverage'] if flag else 0. for flag,r in zip(events,rows[role])])))
            summary.append(row)
            paired.append(dict(condition=condition,method=method,alpha=alpha,lost=lost,gained=gained,
                nominal_p=1. if lost+gained==0 else float(binomtest(gained,lost+gained,.5).pvalue)))
    running=0.
    for rank,i in enumerate(sorted(range(len(paired)),key=lambda i:paired[i]['nominal_p'])):
        running=max(running,min(1.,(len(paired)-rank)*paired[i]['nominal_p']));paired[i]['holm_p']=running
    (a.out/'evaluation.json').write_text(json.dumps(decisions,indent=2))
    (a.out/'summary.json').write_text(json.dumps(dict(summary=summary,retained_historical=retained,family=FAMILY,
        paired=paired,paired_family=27,scope='Corrected numerical wheel, experimental binned pooling/exact path; no default or biological guarantee'),indent=2))


def audit(a,rows,inputs):
    summary=json.loads((a.task/'summary_v1/summary.json').read_text())
    decisions=json.loads((a.task/'summary_v1/evaluation.json').read_text())
    if summary['family']!=99 or len(summary['summary'])!=27 or len(summary['retained_historical'])!=72:
        raise AssertionError('Confidence family changed')
    index={(r['condition'],r['method'],r['alpha'],r['replicate']):r for r in decisions}
    if len(index)!=5400 or len(decisions)!=5400:raise AssertionError('Decision population incomplete')
    count=0;membership={}
    for role,(_,n) in ROLES.items():
        old=read(inputs['old:'+role])
        if len(old)!=n:raise AssertionError('Unpaired old population')
        for r,before in zip(rows[role],old):
            if r['seed_key']!=before['seed_key'] or r['input_sha256']!=before['input_sha256']:
                raise AssertionError('Unpaired observations')
            verify(r['assignments'])
            with np.load(r['assignments']['path']) as z:
                if z['raw_R'].shape[0]!=600 or not np.isfinite(z['raw_R']).all():raise AssertionError('Invalid responsibilities')
                np.testing.assert_array_equal(z['raw_MAP'],z['raw_R'].argmax(axis=1))
                sizes=np.bincount(z['raw_MAP'],minlength=z['raw_R'].shape[1])
                eligible=np.flatnonzero(sizes[1:]>=30)+1
                np.testing.assert_array_equal(eligible,z['eligible_tumour_groups'])
                if len(eligible)!=r['raw_candidate_count']:raise AssertionError('Eligible support differs')
                verify(before['assignments'])
                with np.load(before['assignments']['path']) as old_z:
                    for name in ('simulated_CNA_labels','latent_state'):np.testing.assert_array_equal(z[name],old_z[name])
                if role=='eval_two':
                    truth=z['simulated_CNA_labels'];label=z['raw_MAP']
                    table=[[int(np.count_nonzero((truth==t)&(label==g))) for g in eligible]+[0,0] for t in (0,1)]
                    # Exhaustive two-class injective matching, independent of
                    # the summary's Hungarian implementation. Dummy groups
                    # allow no assignment without credit for extra groups.
                    correct=max(table[0][i]+table[1][j] for i in range(len(table[0])) for j in range(len(table[1])) if i!=j)
                    membership[r['replicate']]=dict(accuracy=correct/600,coverage=int(np.isin(label,eligible).sum())/600,
                                                   groups=len(eligible))
            pair=[];errors=[]
            for d in r['final_pair_diagnostics']:
                values=[abs(x['z']) for x in d['segments']];value=math.fsum(values);pair.append(value)
                eps=np.finfo(float).eps;nterm=len(values)
                errors.append((nterm*eps/(1-nterm*eps)+2*eps)*value)
            if abs(max(pair,default=0.)-r['score'])>max(errors,default=0.):raise AssertionError('Independent final statistic differs')
            if role.startswith('eval_'):
                for method,bank_name,alpha in RULES:
                    bank=None if bank_name is None else [x['score'] for x in rows[bank_name]]
                    prob=None if bank is None else (1+sum(value>=r['score'] for value in bank))/200
                    expected=r['raw_candidate_count']>1 and (bank is None or r['score']>0 and prob<=alpha)
                    saved=index[(r['condition'],method,alpha,r['replicate'])]
                    if saved['event']!=expected or saved['rank_p']!=prob:raise AssertionError('Independent decision mismatch')
            count+=1
    bounds=[];tail=.05/99
    old=json.loads(inputs['old_summary'].read_text())
    old_rows=old['retained_historical']+old['summary']
    for i,r in enumerate(summary['summary']+summary['retained_historical']):
        if i<27:
            population=[x for x in decisions if (x['condition'],x['method'],x['alpha'])==(r['condition'],r['method'],r['alpha'])]
            k,n=sum(x['event'] for x in population),len(population)
        else:k,n=old_rows[i-27]['events'],old_rows[i-27]['n']
        if (k,n)!=(r['events'],r['n']):raise AssertionError('Endpoint count mismatch')
        if i<27 and r['condition']=='two_clones':
            accepted={x['replicate']:x['event'] for x in population}
            if sum(accepted[j] and membership[j]['groups']==2 for j in accepted)!=r['exactly_two_and_accepted']:
                raise AssertionError('Exactly-two power mismatch')
            for metric,name in [('accuracy','mean_correct_fraction_all_cells_if_gate_applied'),
                                ('coverage','mean_coverage_if_gate_applied')]:
                values=[membership[j][metric] if accepted[j] else 0. for j in accepted]
                value=math.fsum(values)/200;eps_mean=np.finfo(float).eps
                allowance=(200*eps_mean/(1-200*eps_mean)+2*eps_mean)*max(value,abs(r[name]))
                if abs(value-r[name])>allowance:raise AssertionError('Independent membership/abstention mismatch')
        eps=np.finfo(float).eps
        if r['direction']=='upper':
            reference=1. if k==n else brentq(lambda p:binom.cdf(k,n,p)-tail,0.,1.,xtol=eps,rtol=4*eps)
        else:reference=0. if k==0 else brentq(lambda p:binom.sf(k-1,n,p)-tail,0.,1.,xtol=eps,rtol=4*eps)
        if abs(reference-r['simultaneous_bound'])>64*eps:raise AssertionError('Independent binomial root differs')
        expected_flag=None if r['direction']=='lower' else reference<=.05
        if r['error_criterion_met']!=expected_flag:raise AssertionError('Bound decision uses wrong confidence family')
        if i>=27 and (r['historical_bound']['family']!=old_rows[i-27]['family'] or
                      r['historical_bound']['simultaneous_bound']!=old_rows[i-27]['simultaneous_bound']):
            raise AssertionError('Historical confidence bound not preserved')
        bounds.append(dict(events=k,n=n,reference=reference))
    prior=json.loads(inputs['old_strict_decisions'].read_text())
    prior += [dict(r,alpha=None if r['method']=='raw' else .05) for r in json.loads(inputs['old_decisions'].read_text())]
    old_index={(r['condition'],r['method'],r['alpha'],r['replicate']):r['event'] for r in prior}
    if len(summary['paired'])!=27:raise AssertionError('Paired family changed')
    for r in summary['paired']:
        keys=[key for key in index if key[:3]==(r['condition'],r['method'],r['alpha'])]
        lost=sum(old_index[key] and not index[key]['event'] for key in keys)
        gained=sum(index[key]['event'] and not old_index[key] for key in keys)
        if len(keys)!=200 or (lost,gained)!=(r['lost'],r['gained']):raise AssertionError('Paired discordance differs')
    snapshot=a.task/'candidate_v1';baseline=json.loads((a.task/'baseline/source_manifest.json').read_text())
    for name,digest in baseline.items():
        if hashlib.sha256((a.task/'baseline'/name).read_bytes()).hexdigest()!=digest:raise AssertionError('Baseline changed')
        if name.startswith(('python/','tests/','crates/')) and hashlib.sha256((snapshot/name).read_bytes()).hexdigest()!=digest:
            raise AssertionError('Executed core changed')
    artifact_path=a.task.parents[1]/'dev/precision_wheel_20261108/validation_22968997/installed_artifact.json'
    artifact=json.loads(artifact_path.read_text());verify(artifact['wheel']);verify(artifact['native'])
    jobs=json.loads((a.task/'jobs.json').read_text());job_ids=[]
    if len(jobs)!=20:raise AssertionError('Wrong execution population')
    for job in jobs:
        identifier=job['result']['output'].strip().split()[-1]
        execution=json.loads((a.task/f'analysis_{identifier}/execution.json').read_text())
        if Path(execution['source']).resolve()!=snapshot.resolve():raise AssertionError('Wrong executed source')
        if execution['wheel']['sha256']!=artifact['wheel']['sha256'] or execution['native_sha256']!=artifact['native']['sha256']:
            raise AssertionError('Wrong executed numerical wheel')
        job_ids.append(identifier)
    (a.out/'report.json').write_text(json.dumps(dict(paired_inputs_and_assignments=count,decisions=5400,
        independent_counts_and_bounds=bounds,paired_discordances=27,baseline_files_preserved=len(baseline),
        independent_membership_power_endpoints=9,corrected_wheel_jobs_verified=job_ids,
        scope='Known-reference simulation audit; not biological validation or a universal clone-membership guarantee'),indent=2))


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--old-task',type=Path,required=True);p.add_argument('--old-tradeoff',type=Path,required=True)
    p.add_argument('--mode',choices=['summary','audit'],required=True);a=p.parse_args()
    rows,inputs=collect(a)
    if a.mode=='audit':
        inputs.update(new_summary=a.task/'summary_v1/summary.json',new_evaluation=a.task/'summary_v1/evaluation.json',
                      job_registry=a.task/'jobs.json',validated_artifact=a.task.parents[1]/'dev/precision_wheel_20261108/validation_22968997/installed_artifact.json')
    a.out.mkdir(exist_ok=False)
    rec=start_run(a.out,inputs,dict(family=99,mode=a.mode),dict(scope='Fixed paired development replay; no parameter selection'),mode='precision_null_'+a.mode)
    try:
        for name,path in inputs.items():
            if name.endswith(':manifest'):
                m=json.loads(path.read_text())
                for record in m['input']['files'].values():verify(record)
                verify(m['output']['files']['records.jsonl'])
        (summaries if a.mode=='summary' else audit)(a,rows,inputs)
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
