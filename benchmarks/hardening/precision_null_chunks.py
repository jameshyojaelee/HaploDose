"""Complete fixed-design precision replay in bounded independent specimen chunks."""
import argparse
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
from haplodose.clones import CloneParams
from haplodose.provenance import start_run,finish_run,file_record


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--role',choices=['cal_nb','cal_state','eval_nb','eval_state','eval_two'],required=True)
    p.add_argument('--start',type=int,required=True);p.add_argument('--stop',type=int,required=True)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--rank-source',type=Path,required=True)
    p.add_argument('--frozen-shards',type=Path,required=True);a=p.parse_args()
    runner=load('unchanged_runner',a.source/'exact_stage_null.py')
    rank=load('unchanged_rank',a.rank_source/'global_rank_gate.py')
    hook=load('unchanged_confirmation',a.source/'shared_membership_stages.py')
    pool=load('unchanged_pooling',Path(__file__).with_name('pooled_initialization.py'))
    gen=rank.load_generator(a.rank_source/'global_rank_generator.py')
    design=json.loads((a.rank_source/'prespecified_design.json').read_text());rank.validate_design(design)
    condition,role,n=runner.ROLES[a.role]
    if a.start not in (0,50,100,150) or a.stop!=min(a.start+50,n):raise ValueError('Unplanned chunk range')
    previous=a.frozen_shards/'exact_matched'/a.role/'records.jsonl'
    originals=[json.loads(line) for line in previous.read_text().splitlines()]
    if len(originals)!=n:raise ValueError('Frozen comparison population changed')
    inputs=dict(script=__file__,pooling=Path(__file__).with_name('pooled_initialization.py'),previous=previous)
    for directory,names in ((a.source,['exact_stage_null.py','shared_membership_stages.py','shared_membership_hmm.py','membership_weighted_stages.py']),
                            (a.rank_source,['global_rank_gate.py','global_rank_generator.py','prespecified_design.json'])):
        for name in names:inputs[name]=directory/name
    params=CloneParams(n_programmes=0,spatial_neighbors=0,seed=design['master_seed'],split_mode='thin')
    a.out.mkdir(parents=True,exist_ok=False)
    rec=start_run(a.out,inputs,dict(params=asdict(params),role=a.role,start=a.start,stop=a.stop,
        design=design,reset_every=20,confidence_family=99,original_alpha=.05,additional_alphas=[.005,.01,.025]),
        dict(scope='Same frozen binned simulation and research pooling/exact path; corrected numerical wheel; no biological guarantee'),
        mode='precision_null_chunk')
    trace=[]
    try:
        with pool.pooled_initialization(trace):
            for i in range(a.start,a.stop):
                first=len(trace)
                row=runner.run_one(rank,gen,hook,design,params,condition,role,i,a.out,'exact_matched')
                if row['seed_key']!=originals[i]['seed_key'] or row['input_sha256']!=originals[i]['input_sha256']:
                    raise AssertionError('Frozen observations changed')
                row['initialization_trace']=trace[first:]
                rank.append(a.out/'records.jsonl',row)
                print(json.dumps(dict(role=a.role,start=a.start,completed=i-a.start+1,total=a.stop-a.start)),flush=True)
        for record in rec['input']['files'].values():
            if file_record(record['path'])['sha256']!=record['sha256']:raise AssertionError('Source/input changed')
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
