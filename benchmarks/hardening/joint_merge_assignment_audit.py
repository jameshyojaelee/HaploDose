"""Recount unchanged criteria and verify the merge assignment intervention."""
import argparse
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record


def verify(record):
    if file_record(record['path'])['sha256']!=record['sha256']:raise AssertionError('Recorded input/output changed')


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--control',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
    helper=Path(__file__).with_name('joint_unordered_audit.py')
    spec=importlib.util.spec_from_file_location('fixed_recount',helper);scorer=importlib.util.module_from_spec(spec);spec.loader.exec_module(scorer)
    candidate=a.task/'candidate_run_v1';inputs=dict(script=Path(__file__),scorer=helper,design=a.task/'DESIGN.md',
        candidate=candidate/'run_manifest.json',control=a.control/'run_manifest.json',assignments=candidate/'assignment_trace.json')
    for seed in (20261031,20261107,20261108,20261109,20261110):
        inputs[f'{seed}:trace']=candidate/f'panel/seed_{seed}/trace.json'
        inputs[f'{seed}:summary']=candidate/f'panel/seed_{seed}/summary.json'
    rec=start_run(a.out,inputs,dict(criteria='original six unchanged'),dict(scope='Selected five-case stage ablation, no calibration claim'),
        mode='joint_merge_assignment_audit')
    try:
        for name in ('candidate','control'):
            m=json.loads(inputs[name].read_text())
            if m['status']!='complete':raise AssertionError('Incomplete inference')
            for r in m['input']['files'].values():verify(r)
            for r in m['output']['files'].values():verify(r)
        cases=[];intercepts=json.loads(inputs['assignments'].read_text());merge_pairs=0
        for item in intercepts:
            full=np.asarray(item['full_responsibility']);prop=np.asarray(item['counterfactual_proposal_responsibility'])
            if item['source']!='joint_full' or full.shape!=prop.shape or len(full)!=len(item['specimen_rows']):raise AssertionError('Invalid intercepted membership')
            if float(np.max(np.abs(full-prop)))!=item['maximum_probability_change']:raise AssertionError('Wrong recorded probability change')
        for seed in (20261031,20261107,20261108,20261109,20261110):
            result=dict(seed=seed)
            for kind,directory in [('control',a.control/'panel'),('candidate',candidate/'panel')]:
                case=directory/f'seed_{seed}';fit=case/'fit';m=json.loads((fit/'run_manifest.json').read_text())
                if m['status']!='complete':raise AssertionError('Incomplete fit')
                for r in m['input']['files'].values():verify(r)
                for r in m['output']['files'].values():verify(r)
                cells=pd.read_csv(fit/'cells.tsv.gz',sep='\t')
                truth=pd.read_csv(a.task.parents[1]/f'dev/recovery_20261107/recovery_panel_v1/seed_{seed}/data/cells.tsv',sep='\t')
                joined=cells.merge(truth[['cell_id','truth_clone']],on='cell_id',validate='one_to_one')
                if len(joined)!=len(cells):raise AssertionError('Scoring lost cells')
                values,passes=scorer.recovery(joined,pd.read_csv(fit/'differences.tsv',sep='\t'))
                recorded=json.loads((case/'summary.json').read_text())
                for key,value in values.items():
                    if recorded['criteria'][key]['value']!=value or recorded['criteria'][key]['passed']!=passes[key]:raise AssertionError('Independent criterion differs')
                result[kind]=dict(values=values,passes=passes,all_numeric_pass=all(passes.values()))
                qc=json.loads((fit/'qc.json').read_text())
                if qc['confirmation_status']!='joint_not_independent':raise AssertionError('False independence claim')
                if kind=='candidate':
                    if qc['merge_assignment_source']!='joint_full':raise AssertionError('Incorrect merge metadata')
                    for entry in json.loads((case/'trace.json').read_text()):
                        if entry['stage']!='merge_membership':continue
                        if entry['assignment_source']!='joint_full' or 'proposal_responsibility' in entry:raise AssertionError('Mislabelled operative membership')
                        np.testing.assert_array_equal(entry['merge_responsibility'],entry['full_responsibility']);merge_pairs+=len(entry['pairs'])
            cases.append(result)
        if not intercepts or not merge_pairs:raise AssertionError('Intervention not exercised')
        baseline=json.loads((a.task/'baseline/source_manifest.json').read_text());core=0
        for name,digest in baseline.items():
            if file_record(a.task/'baseline'/name)['sha256']!=digest:raise AssertionError('Baseline changed')
            if name.startswith(('python/','tests/','crates/')):
                if file_record(a.task.parents[1]/name)['sha256']!=digest:raise AssertionError('Core changed')
                core+=1
        report=dict(cases=cases,intercepted_calls=len(intercepts),verified_operative_merge_pairs=merge_pairs,
            baseline_files=len(baseline),core_files=core,criteria_weakened=False,
            scope='Joint-data research mode; false-split/power and biological validity not established')
        (a.out/'report.json').write_text(json.dumps(report,indent=2));finish_run(a.out,rec,backend='cpu')
        print(json.dumps(dict(candidate_passes=sum(r['candidate']['all_numeric_pass'] for r in cases),
            control_passes=sum(r['control']['all_numeric_pass'] for r in cases),verified_merge_pairs=merge_pairs)),flush=True)
    except Exception as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
