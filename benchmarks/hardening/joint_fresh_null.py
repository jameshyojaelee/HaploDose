"""Fresh, paired complete-binned-workflow comparison of joint merge assignments."""
import argparse
from contextlib import contextmanager,nullcontext
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import jax
from haplodose.clones import CloneParams
from haplodose.provenance import start_run,finish_run,file_record


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--arm',choices=['joint_proposal','joint_full'],required=True)
    p.add_argument('--role',choices=['cal_nb','cal_state','eval_nb','eval_state','eval_two'],required=True)
    p.add_argument('--start',type=int,required=True);p.add_argument('--stop',type=int,required=True)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--rank-source',type=Path,required=True)
    p.add_argument('--joint-source',type=Path,required=True);a=p.parse_args()
    if jax.default_backend()!='cpu':raise RuntimeError('Expected CPU inside inference')
    runner=load('fixed_runner',a.source/'exact_stage_null.py');rank=load('fixed_rank',a.rank_source/'global_rank_gate.py')
    generator=rank.load_generator(a.rank_source/'global_rank_generator.py')
    design=json.loads((a.rank_source/'prespecified_design.json').read_text());rank.validate_design(design)
    condition,old_role,n=runner.ROLES[a.role];role=old_role+100
    if role not in (110,111,120) or role in {x[1] for x in runner.ROLES.values()}:raise AssertionError('Seed role overlaps old panel')
    if a.start not in (0,50,100,150) or a.stop!=min(a.start+50,n):raise ValueError('Unplanned chunk range')
    joint=load('frozen_joint',a.joint_source/'joint_unordered_stages.py')
    full=load('frozen_full_merge',a.joint_source/'joint_merge_assignment.py')
    pool=load('frozen_pooling',a.joint_source/'pooled_initialization.py')
    audit_path=a.task.parents[1]/'dev/joint_merge_assignment_20261108/audit_v1/report.json'
    audit_manifest=json.loads((audit_path.parent/'run_manifest.json').read_text())
    if audit_manifest['status']!='complete':raise AssertionError('Development audit incomplete')
    if file_record(audit_path)['sha256']!=audit_manifest['output']['files']['report.json']['sha256']:raise AssertionError('Development audit changed')
    audited=json.loads(audit_path.read_text())
    if len(audited['cases'])!=5 or not all(r['candidate']['all_numeric_pass'] for r in audited['cases']):raise AssertionError('Recovery gate failed')
    numerical_path=a.task/'preflight_v2/run_manifest.json';numerical=json.loads(numerical_path.read_text())
    if numerical['status']!='complete':raise AssertionError('Numerical preflight incomplete')
    checks=numerical['output']['files']['checks.json']
    if file_record(checks['path'])['sha256']!=checks['sha256']:raise AssertionError('Preflight results changed')
    if file_record(Path(__file__).with_name('joint_fresh_summary.py'))['sha256']!=numerical['input']['files']['statistics']['sha256']:
        raise AssertionError('Statistics changed since numerical preflight')
    @contextmanager
    def stages(stage,trace):
        if stage!=a.arm:raise ValueError('Wrong scientific arm')
        assignment=[]
        with joint.joint_stages('margin',trace),(full.full_joint_merge(assignment) if a.arm=='joint_full' else nullcontext()):
            yield
        trace.extend(dict(stage='full_merge_assignment',**r) for r in assignment)
    hook=SimpleNamespace(weighted_stages=stages)
    inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md',recovery_audit=audit_path,recovery_manifest=audit_path.parent/'run_manifest.json',
        numerical_preflight=numerical_path,numerical_checks=Path(checks['path']),numerical_note=a.task/'NUMERICAL_NOTE_v2.md')
    for directory,names in ((a.source,['exact_stage_null.py']),
            (a.rank_source,['global_rank_gate.py','global_rank_generator.py','prespecified_design.json']),
            (a.joint_source,['joint_unordered_stages.py','joint_merge_assignment.py','shared_membership_stages.py',
                'shared_membership_hmm.py','label_invariant_membership.py','membership_weighted_stages.py','pooled_initialization.py'])):
        for name in names:inputs[name]=directory/name
    params=CloneParams(n_programmes=0,spatial_neighbors=0,seed=design['master_seed'],split_mode='thin')
    a.out.mkdir(parents=True,exist_ok=False)
    rec=start_run(a.out,inputs,dict(params=asdict(params),arm=a.arm,role=a.role,seed_role=role,start=a.start,stop=a.stop,
        generator_design={k:design[k] for k in ('master_seed','cells_per_specimen','bins_per_specimen','generator_sha256')},
        reset_every=20,confidence_family=153,primary_alpha=.05,secondary_alphas=[.005,.01,.025]),
        dict(scope='Fresh conditional known-reference binned workflow; joint-data inference, no biological claim'),mode='joint_fresh_chunk')
    initialization=[]
    try:
        with pool.pooled_initialization(initialization):
            for i in range(a.start,a.stop):
                first=len(initialization)
                row=runner.run_one(rank,generator,hook,design,params,condition,role,i,a.out,a.arm)
                if row['seed_key']!=[20261107,role,rank.CONDITIONS.index(condition),i]:raise AssertionError('Fresh seed changed')
                row['initialization_trace']=initialization[first:];row['arm']=a.arm
                rank.append(a.out/'records.jsonl',row)
                print(json.dumps(dict(arm=a.arm,role=a.role,start=a.start,completed=i-a.start+1,total=a.stop-a.start)),flush=True)
        for r in rec['input']['files'].values():
            if file_record(r['path'])['sha256']!=r['sha256']:raise AssertionError('Source or evidence changed')
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
