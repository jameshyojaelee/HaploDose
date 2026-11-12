"""Five locked cases under explicitly joint-data research inference."""
import argparse
from contextlib import contextmanager
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from haplodose.clones import CloneParams
from haplodose.pipeline import run
from haplodose import pipeline
from haplodose.provenance import start_run, finish_run, file_record


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


@contextmanager
def joint_qc():
    """Label experimental QC before pipeline finish_run hashes it."""
    original=pipeline.write_json
    def write(value,path):
        if Path(path).name=='qc.json':
            value.update(inference_procedure='research_joint_unordered_short_events',confirmation_status='joint_not_independent',
                split_mode_interpretation='Count partition used for proposal search; short-event likelihood reuses full counts',
                original_thinning_assertion_status='not_applicable_to_joint_procedure; original failure retained')
        return original(value,path)
    pipeline.write_json=write
    try:yield
    finally:pipeline.write_json=original


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True);p.add_argument('--panel-root',type=Path,required=True)
    p.add_argument('--assertion-source',type=Path,required=True);p.add_argument('--merge',choices=['hard','margin'],required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);root=Path(__file__).parent
    module=load('original_contract',a.assertion_source);score=load('scoring',root/'recovery_baseline.py')
    hook=load('joint_hook',root/'joint_unordered_stages.py');pool=load('pool_hook',root/'pooled_initialization.py')
    probe=load('label_checks',root/'label_invariance_probe.py')
    model=load('model',root/'shared_membership_hmm.py');unordered=load('unordered',root/'label_invariant_membership.py')
    params=CloneParams(split_mode='thin');frozen=json.loads((a.panel_root/'manifest.json').read_text())
    if tuple(frozen['simulation_seeds'])!=score.SEEDS:raise ValueError('Original seed family changed')
    inputs={path.name:path for path in (root/name for name in (
        'joint_unordered_panel.py','joint_unordered_stages.py','membership_weighted_stages.py','pooled_initialization.py',
        'shared_membership_stages.py','shared_membership_hmm.py','label_invariant_membership.py','label_invariance_probe.py',
        'recovery_baseline.py'))}
    inputs.update(contract=a.assertion_source,panel_manifest=a.panel_root/'manifest.json',design=a.task/'DESIGN.md')
    for entry in frozen['cases']:
        for item in entry['input_files']+[entry['model_metadata']]:
            if file_record(item['path'])['sha256']!=item['sha256']:raise ValueError('Frozen input changed')
            inputs[f"{entry['seed']}:{Path(item['path']).name}"]=Path(item['path'])
    rec=start_run(a.out,inputs,asdict(params),dict(procedure='research_joint_unordered_short_events',merge=a.merge,
        confirmation_status='joint_not_independent',scope='Five exposed development draws, no biological calibration'),
        mode='joint_unordered_recovery')
    try:
        (a.out/'checks.json').write_text(json.dumps(probe.checks(model,unordered),indent=2))
        results=[]
        for seed in score.SEEDS:
            src=a.panel_root/f'seed_{seed}';out=a.out/f'seed_{seed}';out.mkdir();start=time.monotonic()
            metadata=pd.read_csv(src/'model_cells.tsv',sep='\t')
            if any(c.startswith('truth_') for c in metadata):raise ValueError('Truth in inference')
            trace=[];initial=[]
            try:
                with pool.pooled_initialization(initial),hook.joint_stages(a.merge,trace),joint_qc():
                    run(src/'data/counts.h5',src/'model_cells.tsv',src/'data/positions.tsv',out/'fit',chunk=97,
                        params=params,**module.INPUT_CONTRACT)
            finally:
                (out/'trace.json').write_text(json.dumps(trace,indent=2))
                (out/'initialization.json').write_text(json.dumps(initial,indent=2))
            qcpath=out/'fit/qc.json';qc=json.loads(qcpath.read_text())
            if qc['confirmation_status']!='joint_not_independent':raise AssertionError('Unlabelled joint inference')
            fit_manifest=json.loads((out/'fit/run_manifest.json').read_text())
            if fit_manifest['output']['files']['qc.json']['sha256']!=file_record(qcpath)['sha256']:
                raise AssertionError('QC changed after provenance hashing')
            truth=pd.read_csv(src/'data/cells.tsv',sep='\t');fitted=pd.read_csv(out/'fit/cells.tsv.gz',sep='\t')
            joined=fitted.merge(truth[['cell_id']+[c for c in truth if c.startswith('truth_')]],
                on='cell_id',validate='one_to_one')
            joined.to_csv(out/'scored_cells.tsv',sep='\t',index=False)
            checks=score.criteria(joined,pd.read_csv(out/'fit/differences.tsv',sep='\t'))
            row=dict(seed=seed,criteria=checks,all_numeric_criteria_pass=all(x['passed'] for x in checks.values()),
                original_thinning_assertion_status='not_applicable_joint_procedure',seconds=time.monotonic()-start,
                s3_groups=int(joined.loc[joined.specimen.eq('S3') & joined.malignant_call,'clone'].nunique()))
            results.append(row);(out/'summary.json').write_text(json.dumps(row,indent=2));print(json.dumps(row),flush=True)
        (a.out/'summary.json').write_text(json.dumps(dict(merge=a.merge,cases=results,
            independent_confirmation=False,complete_recovery_screen=all(x['all_numeric_criteria_pass'] for x in results)),indent=2))
        baseline=json.loads((a.task/'baseline/source_manifest.json').read_text());core=0
        for name,digest in baseline.items():
            if file_record(a.task/'baseline'/name)['sha256']!=digest:raise AssertionError('Baseline changed')
            if name.startswith(('python/','tests/','crates/')):
                if file_record(a.task.parents[1]/name)['sha256']!=digest:raise AssertionError('Core changed')
                core+=1
        for r in rec['input']['files'].values():
            if file_record(r['path'])['sha256']!=r['sha256']:raise AssertionError('Input changed')
        (a.out/'preservation.json').write_text(json.dumps(dict(baseline_files=len(baseline),core_files=core,
            inputs=len(rec['input']['files'])),indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
