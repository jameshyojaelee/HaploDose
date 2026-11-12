"""Recount activated-workflow endpoints and independently invert binomial tails."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.stats import binom


def read(path):return [json.loads(line) for line in path.read_text().splitlines()]


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    summary=json.loads((a.task/'null_summary_v1/summary.json').read_text())
    evaluation=json.loads((a.task/'null_summary_v1/evaluation.json').read_text());root=a.task.parents[1]
    historical_evaluation=json.loads((root/'dev/confirmation_20261107/null_summary_v1/evaluation.json').read_text())
    checks=[];input_checks=0;mechanism={}
    for stage in ('merge_control','exact_matched'):
        for role,n in [('cal_nb',199),('cal_state',199),('eval_nb',200),('eval_state',200),('eval_two',200)]:
            path=a.task/'null_shards_v1'/stage/role
            manifest=json.loads((path/'run_manifest.json').read_text())
            if manifest['status']!='complete':raise ValueError('Incomplete shard')
            for record in manifest['input']['files'].values():
                if hashlib.sha256(Path(record['path']).read_bytes()).hexdigest()!=record['sha256']:
                    raise ValueError('Executed source changed')
            current=read(path/'records.jsonl');old=read(root/'dev/confirmation_20261107/null_shards_v1/merge'/role/'records.jsonl')
            if len(current)!=n or len(old)!=n:raise ValueError('Missing original or current replicates')
            for before,after in zip(old,current):
                if before['seed_key']!=after['seed_key']:raise ValueError('Seed changed')
                for key,h in before['input_sha256'].items():
                    if key!='reset' and after['input_sha256'][key]!=h:raise ValueError('Counts/nuisance generator changed')
                if before['input_sha256']['reset']==after['input_sha256']['reset']:raise ValueError('Resolution challenge not active')
                input_checks+=1
            bank_nb=[r['score'] for r in read(a.task/'null_shards_v1'/stage/'cal_nb/records.jsonl')]
            bank_state=[r['score'] for r in read(a.task/'null_shards_v1'/stage/'cal_state/records.jsonl')]
            for row in current:
                if role.startswith('eval_'):
                    for method,bank in [('raw',None),('NB_calibrated_gate',bank_nb),('state_calibrated_gate',bank_state)]:
                        probability=None if bank is None else (1+sum(value>=row['score'] for value in bank))/200
                        event=bool(row['raw_candidate_count']>1 and (bank is None or row['score']>0 and probability<=.05))
                        saved=[r for r in evaluation if r['stage']==stage and r['condition']==row['condition']
                               and r['method']==method and r['replicate']==row['replicate']]
                        if len(saved)!=1 or saved[0]['event']!=event or saved[0]['rank_p']!=probability:
                            raise ValueError('Independent decision mismatch')
            calls=sum(len(e.get('calls',[])) for r in current for e in r['stage_trace'] if e['stage']=='exact_confirmation')
            mechanism[stage+':'+role]=calls
    if not sum(n for name,n in mechanism.items() if name.startswith('exact_matched:')):
        raise ValueError('Candidate mechanism was never exercised')
    alpha=.05/45
    count_checks=[]
    for historical,group in [(False,summary['summary']),(True,summary['retained_historical'])]:
        population=historical_evaluation if historical else evaluation
        for row in group:
            selected=[r for r in population if all(r[key]==row[key] for key in ('stage','condition','method'))]
            count=sum(r['event'] for r in selected)
            if len(selected)!=row['n'] or count!=row['events']:
                raise ValueError('Summary endpoint count differs from independently checked decisions')
            count_checks.append(dict(historical=historical,stage=row['stage'],condition=row['condition'],
                                     method=row['method'],events=count,n=len(selected)))
    for row in summary['paired']:
        selected=[[r for r in evaluation if r['stage']==stage and r['condition']==row['condition'] and r['method']==row['method']]
                  for stage in ('merge_control','exact_matched')]
        before,after=[{r['replicate']:r['event'] for r in rows} for rows in selected]
        if len(before)!=200 or before.keys()!=after.keys():raise ValueError('Unpaired specimen decisions')
        lost=sum(before[i] and not after[i] for i in before);gained=sum(after[i] and not before[i] for i in before)
        if (lost,gained)!=(row['lost'],row['gained']):raise ValueError('Discordance count mismatch')
    for row in summary['summary']+summary['retained_historical']:
        k,n=row['events'],row['n']
        if row['direction']=='upper':
            ref=1. if k==n else brentq(lambda x:binom.cdf(k,n,x)-alpha,0.,1.,xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps)
        else:
            ref=0. if k==0 else brentq(lambda x:binom.sf(k-1,n,x)-alpha,0.,1.,xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps)
        discrepancy=abs(ref-row['simultaneous_bound'])
        if discrepancy>64*np.finfo(float).eps:raise ArithmeticError('Independent bound disagreement')
        checks.append(dict(stage=row['stage'],condition=row['condition'],method=row['method'],events=k,n=n,
                           bound=ref,absolute_discrepancy=discrepancy))
    report=dict(preserved_generator_inputs=input_checks,independent_decisions=3600,independent_bounds=checks,
                independently_recounted_endpoints=count_checks,paired_discordances_recounted=len(summary['paired']),
                mechanism_calls=mechanism,scope='Conditional simulated workflows; no biological replication or calibration claim')
    (a.out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)


if __name__=='__main__':main()
