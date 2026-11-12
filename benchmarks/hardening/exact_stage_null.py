"""Fixed resolution challenge that exercises exact confirmation in full bin search.

Counts and nuisance laws are the previous frozen generator, with chromosome
blocks of 20 bins (ten-bin half-chromosome proposals). This is a separate
development challenge, not a replacement for the original 60-bin-block results.
"""
import argparse
from collections import Counter
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
import numpy as np
from scipy.stats import beta,binomtest
from sklearn.metrics import adjusted_rand_score
from haplodose.clones import CloneParams,Rows,fit_specimen
from haplodose.provenance import start_run,finish_run,file_record

ROLES={'cal_nb':('matched_nb',10,199),'cal_state':('correlated_state',11,199),
       'eval_nb':('matched_nb',20,200),'eval_state':('correlated_state',20,200),'eval_two':('two_clones',20,200)}
STAGES=('merge_control','exact_matched')
METHODS=('raw','NB_calibrated_gate','state_calibrated_gate')
FAMILY=45  # original27 + new18; do not discard earlier failures


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def run_one(m,gen,hook,design,params,condition,role,i,out,stage):
    key=[design['master_seed'],role,m.CONDITIONS.index(condition),i]
    rng=np.random.default_rng(np.random.SeedSequence(key))
    Y,base,phi,truth,latent=gen.generate(rng,condition,design['cells_per_specimen'],design['bins_per_specimen'],return_latent=True)
    reset=np.zeros(Y.shape[1],bool);reset[::20]=True
    hashes={name:m.array_hash(value) for name,value in dict(counts=Y,reference=base,phi=phi,truth=truth,reset=reset,
             state_effect=latent['state_effect'],mean=latent['mean']).items()}
    trace=[]
    with hook.weighted_stages(stage,trace):
        result=fit_specimen(Rows(Y,np.arange(len(Y))),Y.sum(axis=1,dtype=np.float64),base[:,None],reset,phi,params)
        score,labels,eligible,differences=m.final_difference_score(result,params.min_cells)
    path=out/f'{condition}_{role}_{i:03d}.npz'
    np.savez_compressed(path,simulated_CNA_labels=truth,raw_MAP=labels,raw_R=result['R'],
        latent_state=latent['state_effect'],eligible_tumour_groups=eligible)
    return dict(condition=condition,role=role,replicate=i,seed_key=key,input_sha256=hashes,score=score,
        raw_candidate_count=len(eligible),all_MAP_group_sizes=dict(Counter(labels.tolist())),
        ARI_against_simulated_CNA_labels=float(adjusted_rand_score(truth,labels)),
        reference_like_fraction=float(np.mean(labels==0)),counted_tumour_fraction=float(np.isin(labels,eligible).mean()),
        assignments=file_record(path),final_pair_diagnostics=differences,stage_trace=trace)


def summarize(a,m):
    a.out.mkdir(exist_ok=False);all_rows={};coverage=[]
    for stage in STAGES:
        all_rows[stage]={}
        for role,(condition,role_id,n) in ROLES.items():
            directory=a.shard_root/stage/role
            rows=[json.loads(x) for x in (directory/'records.jsonl').read_text().splitlines()]
            manifest=json.loads((directory/'run_manifest.json').read_text())
            if len(rows)!=n or manifest['status']!='complete':raise ValueError('Incomplete scientific shard')
            count=sum(len(e.get('calls',[])) for r in rows for e in r['stage_trace'] if e['stage']=='exact_confirmation')
            coverage.append(dict(stage=stage,role=role,exact_likelihood_calls=count))
            for i,row in enumerate(rows):
                if row['seed_key']!=[20261107,role_id,m.CONDITIONS.index(condition),i]:raise ValueError('Seed/order changed')
            all_rows[stage][role]=rows
    evaluations=[];summary=[];alpha=.05/FAMILY
    def bound(k,n,power):
        return (0. if not k else float(beta.ppf(alpha,k,n-k+1))) if power else (1. if k==n else float(beta.ppf(1-alpha,k+1,n-k)))
    for stage,roles in all_rows.items():
        for role in ('eval_nb','eval_state','eval_two'):
            condition=ROLES[role][0]
            for method in METHODS:
                bank=None if method=='raw' else [r['score'] for r in roles['cal_nb' if method=='NB_calibrated_gate' else 'cal_state']]
                events=[]
                for row in roles[role]:
                    probability=None if bank is None else m.rank_p(bank,row['score'])
                    event=bool(row['raw_candidate_count']>1 and (bank is None or (row['score']>0 and probability<=.05)))
                    events.append(event);evaluations.append(dict(stage=stage,condition=condition,method=method,
                        replicate=row['replicate'],event=event,rank_p=probability,score=row['score'],raw_groups=row['raw_candidate_count']))
                k,n=sum(events),len(events);power=condition=='two_clones'
                summary.append(dict(stage=stage,condition=condition,method=method,events=k,n=n,
                    direction='lower' if power else 'upper',simultaneous_bound=bound(k,n,power),
                    raw_ARI_median=float(np.median([r['ARI_against_simulated_CNA_labels'] for r in roles[role]])),
                    raw_exactly_two=sum(r['raw_candidate_count']==2 for r in roles[role])))
    paired=[]
    for role in ROLES:
        for arow,brow in zip(all_rows[STAGES[0]][role],all_rows[STAGES[1]][role]):
            if arow['input_sha256']!=brow['input_sha256']:raise ValueError('Unpaired inputs')
    for condition in m.CONDITIONS:
        for method in METHODS:
            x,y=[[r['event'] for r in evaluations if r['stage']==s and r['condition']==condition and r['method']==method] for s in STAGES]
            lost=sum(a and not b for a,b in zip(x,y));gained=sum(b and not a for a,b in zip(x,y))
            paired.append(dict(condition=condition,method=method,lost=lost,gained=gained,
                nominal_p=1. if lost+gained==0 else float(binomtest(gained,lost+gained,.5).pvalue)))
    largest=0.
    for rank,i in enumerate(sorted(range(len(paired)),key=lambda i:paired[i]['nominal_p'])):
        largest=max(largest,min(1.,(len(paired)-rank)*paired[i]['nominal_p']));paired[i]['holm_p']=largest
    old=json.loads((a.previous/'summary.json').read_text())['summary']
    if len(old)!=27:raise ValueError('Original confidence family incomplete')
    retained=[dict(row,simultaneous_bound=bound(row['events'],row['n'],row['direction']=='lower'),family=FAMILY) for row in old]
    (a.out/'evaluation.json').write_text(json.dumps(evaluations,indent=2))
    (a.out/'summary.json').write_text(json.dumps(dict(summary=summary,retained_historical=retained,paired=paired,
        mechanism_coverage=coverage,confidence_family=FAMILY,tail_alpha=alpha,
        interpretation='Known-reference expression-only full binned workflow; new resolution challenge; not biological validation or membership certification'),indent=2))
    print(json.dumps(summary),flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--rank-source',type=Path,required=True)
    p.add_argument('--shard-root',type=Path,required=True);p.add_argument('--stage',choices=STAGES,default='exact_matched')
    p.add_argument('--role',choices=list(ROLES)+['summary'],required=True);p.add_argument('--previous',type=Path)
    a=p.parse_args();m=load('fixed_rank',a.rank_source/'global_rank_gate.py')
    design=json.loads((a.rank_source/'prespecified_design.json').read_text());m.validate_design(design)
    if a.role=='summary':summarize(a,m);return
    gen=m.load_generator(a.rank_source/'global_rank_generator.py');root=Path(__file__).parent
    hook=load('exact_stages',root/'shared_membership_stages.py')
    params=CloneParams(n_programmes=0,spatial_neighbors=0,seed=design['master_seed'],split_mode='thin')
    a.out.mkdir(parents=True,exist_ok=False)
    inputs={name:root/name for name in ('exact_stage_null.py','shared_membership_stages.py','shared_membership_hmm.py','membership_weighted_stages.py')}
    inputs.update({name:a.rank_source/name for name in ('global_rank_gate.py','global_rank_generator.py','prespecified_design.json')})
    rec=start_run(a.out,inputs,dict(params=asdict(params),stage=a.stage,role=a.role,design=design,chromosome_block_bins=20),
        dict(evidence='new resolution diagnostic on existing generator, fixed original counts law',primary_family=FAMILY,
             scope='full bin search, known reference; public preprocessing, gene and allele vetoes absent'),mode='exact_stage_null')
    try:
        condition,role,n=ROLES[a.role]
        for i in range(n):
            row=run_one(m,gen,hook,design,params,condition,role,i,a.out,a.stage)
            m.append(a.out/'records.jsonl',row)
            print(json.dumps(dict(stage=a.stage,role=a.role,completed=i+1,total=n)),flush=True)
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
