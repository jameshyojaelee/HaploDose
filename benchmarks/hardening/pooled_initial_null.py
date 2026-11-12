"""Complete workflow challenge for pooling initial tumour groups.

Fixed old observation law, seed keys, operating rules and replicate counts.
This is development calibration, not biological validation.
"""
import argparse
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path

from scipy import stats
from haplodose.clones import CloneParams
from haplodose.provenance import start_run,finish_run

FAMILY=54  # retained45 plus new9, no failed endpoint removed


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def bound(k,n,power):
    alpha=.05/FAMILY
    return (0. if k==0 else float(stats.beta.ppf(alpha,k,n-k+1))) if power else (
        1. if k==n else float(stats.beta.ppf(1-alpha,k+1,n-k)))


def summarize(a,runner,rank):
    rows={};inputs={};coverage=[]
    for role,(_,_,n) in runner.ROLES.items():
        directory=a.shards/role
        manifest=json.loads((directory/'run_manifest.json').read_text())
        records=[json.loads(line) for line in (directory/'records.jsonl').read_text().splitlines()]
        if manifest['status']!='complete' or len(records)!=n:raise AssertionError('Incomplete shard')
        rows[role]=records;inputs[role]=directory/'records.jsonl';inputs[role+'_manifest']=directory/'run_manifest.json'
        coverage.append(dict(role=role,exact_calls=sum(len(e.get('calls',[])) for r in records
                        for e in r['stage_trace'] if e['stage']=='exact_confirmation')))
    inputs['old_summary']=a.previous/'summary.json';inputs['old_evaluation']=a.previous/'evaluation.json'
    inputs['script']=Path(__file__)
    a.out.mkdir(exist_ok=False)
    rec=start_run(a.out,inputs,dict(family=FAMILY,operating_alpha=.05),
                  dict(scope='Development full binned workflow with known reference'),mode='pooled_initial_null_summary')
    try:
        evaluation=[];summary=[]
        for role in ('eval_nb','eval_state','eval_two'):
            condition=runner.ROLES[role][0]
            for method in runner.METHODS:
                bank=None if method=='raw' else [r['score'] for r in rows['cal_nb' if method=='NB_calibrated_gate' else 'cal_state']]
                events=[]
                for r in rows[role]:
                    p=None if bank is None else rank.rank_p(bank,r['score'])
                    event=bool(r['raw_candidate_count']>1 and (bank is None or (r['score']>0 and p<=.05)))
                    events.append(event)
                    evaluation.append(dict(stage='pooled_initial',condition=condition,method=method,
                                           replicate=r['replicate'],event=event,rank_p=p,score=r['score'],
                                           raw_groups=r['raw_candidate_count']))
                k,n=sum(events),len(events);power=condition=='two_clones'
                summary.append(dict(stage='pooled_initial',condition=condition,method=method,events=k,n=n,
                                    direction='lower' if power else 'upper',simultaneous_bound=bound(k,n,power),
                                    raw_exactly_two=sum(r['raw_candidate_count']==2 for r in rows[role]),
                                    raw_ARI_mean=sum(r['ARI_against_simulated_CNA_labels'] for r in rows[role])/n))
        old=json.loads(inputs['old_summary'].read_text());prior=old['summary']+old['retained_historical']
        if len(prior)!=45:raise AssertionError('Historical family incomplete')
        retained=[dict(r,simultaneous_bound=bound(r['events'],r['n'],r['direction']=='lower'),family=FAMILY) for r in prior]
        old_decisions=json.loads(inputs['old_evaluation'].read_text());paired=[]
        for stage in runner.STAGES:
            for condition in rank.CONDITIONS:
                for method in runner.METHODS:
                    before=[r['event'] for r in old_decisions if r['stage']==stage and r['condition']==condition and r['method']==method]
                    after=[r['event'] for r in evaluation if r['condition']==condition and r['method']==method]
                    if len(before)!=200 or len(after)!=200:raise AssertionError('Paired population incomplete')
                    lost=sum(x and not y for x,y in zip(before,after));gained=sum(y and not x for x,y in zip(before,after))
                    paired.append(dict(comparator=stage,condition=condition,method=method,lost=lost,gained=gained,
                                       nominal_p=1. if not lost+gained else float(stats.binomtest(gained,lost+gained,.5).pvalue)))
        largest=0.
        for rank_index,i in enumerate(sorted(range(len(paired)),key=lambda i:paired[i]['nominal_p'])):
            largest=max(largest,min(1.,(len(paired)-rank_index)*paired[i]['nominal_p']));paired[i]['holm_p']=largest
        (a.out/'evaluation.json').write_text(json.dumps(evaluation,indent=2))
        (a.out/'summary.json').write_text(json.dumps(dict(summary=summary,retained_historical=retained,
            confidence_family=FAMILY,tail_alpha=.05/FAMILY,paired=paired,paired_family=18,
            mechanism_coverage=coverage,scope='Expression-only binned workflow; no gene/allele/public preprocessing guarantee'),indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--role',choices=['cal_nb','cal_state','eval_nb','eval_state','eval_two','summary'],required=True)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--rank-source',type=Path,required=True)
    p.add_argument('--frozen-shards',type=Path);p.add_argument('--shards',type=Path);p.add_argument('--previous',type=Path)
    a=p.parse_args();runner=load('locked_runner',a.source/'exact_stage_null.py')
    rank=load('locked_rank',a.rank_source/'global_rank_gate.py')
    if a.role=='summary':summarize(a,runner,rank);return
    intervention=load('pooling_intervention',Path(__file__).with_name('pooled_initialization.py'))
    hook=load('locked_exact_hooks',a.source/'shared_membership_stages.py')
    gen=rank.load_generator(a.rank_source/'global_rank_generator.py')
    design=json.loads((a.rank_source/'prespecified_design.json').read_text());rank.validate_design(design)
    params=CloneParams(n_programmes=0,spatial_neighbors=0,seed=design['master_seed'],split_mode='thin')
    a.out.mkdir(parents=True,exist_ok=False)
    inputs={'script':Path(__file__),'intervention':Path(__file__).with_name('pooled_initialization.py')}
    for directory,names in ((a.source,['exact_stage_null.py','shared_membership_stages.py','shared_membership_hmm.py','membership_weighted_stages.py']),
                            (a.rank_source,['global_rank_gate.py','global_rank_generator.py','prespecified_design.json'])):
        for name in names:inputs[name]=directory/name
    inputs['previous_records']=a.frozen_shards/'exact_matched'/a.role/'records.jsonl'
    rec=start_run(a.out,inputs,dict(params=asdict(params),role=a.role,design=design,reset_every=20,
                                  initialization='pool_initial_tumour_groups',confidence_family=FAMILY),
                  dict(scope='Complete binned expression search; no gene or allele observations'),mode='pooled_initial_null')
    trace=[]
    try:
        condition,role_id,n=runner.ROLES[a.role]
        with inputs['previous_records'].open() as old, intervention.pooled_initialization(trace):
            for i in range(n):
                previous=json.loads(next(old));first=len(trace)
                row=runner.run_one(rank,gen,hook,design,params,condition,role_id,i,a.out,'exact_matched')
                if row['input_sha256']!=previous['input_sha256'] or row['seed_key']!=previous['seed_key']:
                    raise AssertionError('Original paired input changed')
                row['initialization_trace']=trace[first:]
                rank.append(a.out/'records.jsonl',row)
                print(json.dumps(dict(role=a.role,completed=i+1,total=n)),flush=True)
            if old.read().strip():raise AssertionError('Original input population changed')
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
