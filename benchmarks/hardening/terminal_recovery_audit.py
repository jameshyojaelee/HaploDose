"""Independent original-criterion recount and terminal-reference parity."""
import argparse
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record


def verify(record):
    if file_record(record['path'])['sha256']!=record['sha256']:raise AssertionError('Recorded artifact changed')


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);root=a.task.parents[1]
    helper=Path(__file__).with_name('joint_unordered_audit.py')
    spec=importlib.util.spec_from_file_location('independent_recount',helper);score=importlib.util.module_from_spec(spec);spec.loader.exec_module(score)
    previous=dict(original=root/'dev/hardware_trace_20261108/recovery_default_v1/panel',
        pooled_exact=root/'dev/hardware_trace_20261108/recovery_pooled_exact_v1/panel',
        joint=root/'dev/joint_unordered_20261108/margin_v2')
    reference=root/'dev/terminal_merge_20261108/probe_v1/panel'
    inputs=dict(script=Path(__file__),scorer=helper,design=a.task/'DESIGN.md')
    for mode in previous:
        inputs[mode]=a.task/f'{mode}_v1/run_manifest.json'
        inputs[mode+':before']=previous[mode]/'summary.json'
        for seed in (20261031,20261107,20261108,20261109,20261110):
            for prefix,directory in [('new',a.task/f'{mode}_v1/panel'),('previous',previous[mode])]:
                case=directory/f'seed_{seed}'
                inputs[f'{mode}:{seed}:{prefix}:summary']=case/'summary.json'
                inputs[f'{mode}:{seed}:{prefix}:manifest']=case/'fit/run_manifest.json'
                if mode=='joint' and prefix=='new':inputs[f'{mode}:{seed}:terminal_reference']=reference/f'seed_{seed}/fit/run_manifest.json'
    rec=start_run(a.out,inputs,dict(modes=list(previous)),dict(scope='Fixed recovery criteria and independent terminal reference, not biological validation'),mode='terminal_recovery_audit')
    try:
        results=[];parity=0;stops=[]
        for mode in previous:
            manifest=json.loads(inputs[mode].read_text())
            if manifest['status']!='complete':raise AssertionError('Incomplete recovery job')
            for r in manifest['input']['files'].values():verify(r)
            panel=a.task/f'{mode}_v1/panel';summary=json.loads((panel/'summary.json').read_text())['cases']
            for row in summary:
                seed=row['seed'];case=panel/f'seed_{seed}';fit=case/'fit'
                fm=json.loads((fit/'run_manifest.json').read_text())
                if fm['status']!='complete':raise AssertionError('Incomplete pipeline run')
                for r in fm['input']['files'].values():verify(r)
                for r in fm['output']['files'].values():verify(r)
                exported=pd.read_csv(fit/'cells.tsv.gz',sep='\t')
                truth=pd.read_csv(root/f'dev/recovery_20261107/recovery_panel_v1/seed_{seed}/data/cells.tsv',sep='\t')
                joined=exported.merge(truth[['cell_id','truth_clone']],on='cell_id',validate='one_to_one')
                if len(joined)!=len(exported):raise AssertionError('Scoring lost cells')
                values,passes=score.recovery(joined,pd.read_csv(fit/'differences.tsv',sep='\t'))
                for key,value in values.items():
                    if row['criteria'][key]['value']!=value or row['criteria'][key]['passed']!=passes[key]:raise AssertionError('Recovery recount differs')
                before=json.loads((previous[mode]/f'seed_{seed}/summary.json').read_text())['criteria']
                results.append(dict(mode=mode,seed=seed,before=before,after=values,passes=passes))
                log=pd.read_csv(fit/'log.tsv',sep='\t')
                for stop in log.loc[log.action=='search_stop'].to_dict('records'):
                    if stop['reason']!='cycle_budget_exhausted' or not stop['terminal_merge_checked']:raise AssertionError('Missing stopping provenance')
                    stops.append(dict(mode=mode,seed=seed,specimen=stop['specimen'],rounds=int(stop['terminal_merge_rounds'])))
                if mode=='joint':
                    target=reference/f'seed_{seed}/fit';tm=json.loads((target/'run_manifest.json').read_text())
                    for r in tm['output']['files'].values():verify(r)
                    np.testing.assert_array_equal(np.load(fit/'posteriors.npy'),np.load(target/'posteriors.npy'))
                    for name in ('cells.tsv.gz','profiles.tsv.gz','differences.tsv','segments.tsv'):
                        pd.testing.assert_frame_equal(pd.read_csv(fit/name,sep='\t'),pd.read_csv(target/name,sep='\t'),check_exact=True)
                    parity+=1
                    if json.loads((fit/'qc.json').read_text())['confirmation_status']!='joint_not_independent':raise AssertionError('Joint mislabelled')
        baseline=json.loads((a.task/'baseline/source_manifest.json').read_text());core=0
        for name,digest in baseline.items():
            if file_record(a.task/'baseline'/name)['sha256']!=digest:raise AssertionError('Baseline changed')
            if name.startswith(('python/','tests/','crates/')):
                if file_record(root/name)['sha256']!=digest:raise AssertionError('Candidate core changed')
                core+=1
        (a.out/'report.json').write_text(json.dumps(dict(cases=results,terminal_reference_exact_cases=parity,
            stopping_records=stops,baseline_files=len(baseline),core_files=core),indent=2))
        finish_run(a.out,rec,backend='cpu');print(json.dumps(dict(cases=len(results),terminal_reference_exact_cases=parity)),flush=True)
    except Exception as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
