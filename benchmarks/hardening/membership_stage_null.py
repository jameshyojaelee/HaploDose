"""Fixed-design complete binned-workflow replay of weighted stage interventions.

The original final-group score includes initial groups and accepted proposals.
Known-reference, expression-only simulation: not biological clone validation.
"""
import argparse
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
import numpy as np
from scipy.stats import beta, binomtest
from haplodose.clones import CloneParams

ROLES={'cal_nb':('matched_nb',10,199),'cal_state':('correlated_state',11,199),
       'eval_nb':('matched_nb',20,200),'eval_state':('correlated_state',20,200),
       'eval_two':('two_clones',20,200)}
STAGES=('hard','merge','both')
METHODS=('raw','NB_calibrated_gate','state_calibrated_gate')

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def read_shards(root,stage,m,params,design):
    records={}
    for name,(condition,role,n) in ROLES.items():
        path=root/stage/name
        meta=json.loads((path/'shard.json').read_text())
        if meta['params']!=json.loads(json.dumps(asdict(params))) or meta['design']!=design or meta['stage']!=stage:
            raise ValueError('Mismatched shard configuration')
        rows=[json.loads(l) for l in (path/'records.jsonl').read_text().splitlines()]
        if len(rows)!=n:raise ValueError('Incomplete shard')
        for i,row in enumerate(rows):
            if row['seed_key']!=[design['master_seed'],role,m.CONDITIONS.index(condition),i]:
                raise ValueError('Wrong specimen identity/order')
        records[name]=rows
    return records

def summarize(a,m,params,design):
    a.out.mkdir(exist_ok=False)
    all_records={s:read_shards(a.shard_root,s,m,params,design) for s in STAGES}
    joint_alpha=.05/(len(STAGES)*len(METHODS)*3)
    summary=[];evaluations=[];paired=[]
    for stage,records in all_records.items():
        for role in ('eval_nb','eval_state','eval_two'):
            condition=ROLES[role][0]
            for method in METHODS:
                events=[]
                bank=None if method=='raw' else [r['score'] for r in records['cal_nb' if method=='NB_calibrated_gate' else 'cal_state']]
                for row in records[role]:
                    p=None if bank is None else m.rank_p(bank,row['score'])
                    event=bool(row['raw_candidate_count']>1 and (bank is None or (row['score']>0 and p<=.05)))
                    events.append(event)
                    evaluations.append(dict(stage=stage,condition=condition,method=method,replicate=row['replicate'],
                                            event=event,rank_p=p,score=row['score'],raw_groups=row['raw_candidate_count']))
                k,n=sum(events),len(events);power=condition=='two_clones'
                bound=(0. if k==0 else float(beta.ppf(joint_alpha,k,n-k+1))) if power else (
                    1. if k==n else float(beta.ppf(1-joint_alpha,k+1,n-k)))
                summary.append(dict(stage=stage,condition=condition,method=method,events=k,n=n,
                                    direction='lower' if power else 'upper',simultaneous_bound=bound,
                                    raw_ARI_median=float(np.median([r['ARI_against_simulated_CNA_labels'] for r in records[role]])),
                                    raw_exactly_two=sum(r['raw_candidate_count']==2 for r in records[role])))
    for stage in STAGES[1:]:
        for role in ROLES:
            for old,new in zip(all_records['hard'][role],all_records[stage][role]):
                if old['input_sha256']!=new['input_sha256']:raise ValueError('Unpaired generated inputs')
        for condition in m.CONDITIONS:
            for method in METHODS:
                rows=[[r['event'] for r in evaluations if r['stage']==s and r['condition']==condition and r['method']==method] for s in ('hard',stage)]
                lost=sum(x and not y for x,y in zip(*rows));gained=sum(y and not x for x,y in zip(*rows))
                p=1. if lost+gained==0 else float(binomtest(gained,lost+gained,.5).pvalue)
                paired.append(dict(stage=stage,condition=condition,method=method,lost=lost,gained=gained,nominal_p=p))
    largest=0.
    for rank,i in enumerate(sorted(range(len(paired)),key=lambda i:paired[i]['nominal_p'])):
        largest=max(largest,min(1.,(len(paired)-rank)*paired[i]['nominal_p']));paired[i]['holm_p']=largest
    (a.out/'evaluation.json').write_text(json.dumps(evaluations,indent=2))
    (a.out/'summary.json').write_text(json.dumps(dict(summary=summary,paired=paired,
        joint_confidence_family=27,tail_alpha=joint_alpha,params=asdict(params),design=design,
        interpretation='Development regression; conditional on realized calibration banks. A global RNA distinction does not validate membership. Known-reference binned workflow excludes public preprocessing, genes and alleles.'),indent=2))
    print(json.dumps(summary),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--rank-source',type=Path,required=True);p.add_argument('--shard-root',type=Path,required=True)
    p.add_argument('--stage',choices=STAGES,default='both');p.add_argument('--role',choices=list(ROLES)+['summary'],required=True)
    a=p.parse_args();m=load('unchanged_rank',a.rank_source/'global_rank_gate.py')
    design=json.loads((a.rank_source/'prespecified_design.json').read_text());m.validate_design(design)
    params=CloneParams(n_programmes=0,spatial_neighbors=0,seed=design['master_seed'],split_mode='thin')
    if a.role=='summary':summarize(a,m,params,design);return
    hook=load('weighted_stages',Path(__file__).with_name('membership_weighted_stages.py'))
    gen=m.load_generator(a.rank_source/'global_rank_generator.py')
    a.out.mkdir(exist_ok=False);condition,role,n=ROLES[a.role]
    (a.out/'shard.json').write_text(json.dumps(dict(stage=a.stage,role=a.role,params=asdict(params),design=design,
        sources=[m.record(__file__),m.record(Path(__file__).with_name('membership_weighted_stages.py'))]+
        [m.record(a.rank_source/f) for f in ('global_rank_gate.py','global_rank_generator.py','prespecified_design.json')]),indent=2))
    for i in range(n):
        trace=[]
        with hook.weighted_stages(a.stage,trace):
            row=m.run_specimen(gen,design,params,condition,role,i,a.out)
        row['stage_trace']=trace
        m.append(a.out/'records.jsonl',row)
        print(json.dumps(dict(stage=a.stage,role=a.role,completed=i+1,total=n)),flush=True)

if __name__=='__main__':main()
