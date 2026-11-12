"""Diagnose unvalidated final splits without changing the installed package."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
from haplodose import pipeline
from haplodose.provenance import start_run,finish_run,file_record


@contextmanager
def terminal_merge(trace):
    original=pipeline.fit_specimen
    def fit(*args,**kwargs):
        result=original(*args,**kwargs);model=result['model']
        if model.alleles is not None or model.graph is not None:raise ValueError('Expression-only nonspatial probe')
        record=dict(specimen_rows=model.rows.rows.tolist(),initial_groups=len(result['doses']),
            before_R=result['R'].tolist(),before_pi=result['pi'].tolist(),checks=[])
        trace.append(record);current=result;limit=len(result['doses'])
        for iteration in range(limit):
            before=len(current['doses']);log_start=len(model.log)
            doses,pi,merged=model.merge_unsupported(current['R'],current['doses'],current['pi'])
            record['checks'].append(dict(iteration=iteration,before_groups=before,after_groups=len(doses),
                merged=bool(merged),actions=model.log[log_start:]))
            if not merged:
                # A second application must be a no-op on returned quantities.
                again_doses,again_pi,again_merged=model.merge_unsupported(current['R'],doses,pi)
                if again_merged:raise AssertionError('Terminal merge check not stable')
                np.testing.assert_array_equal(doses,again_doses)
                np.testing.assert_array_equal(pi,again_pi)
                break
            if len(doses)>=before:raise AssertionError('Merge failed to reduce number of groups')
            doses,pi,R,LL,Yk,Lk,XLk=model.em(doses,pi)
            if len(doses)>=before:raise AssertionError('Terminal refit added groups')
            current=dict(current,doses=np.asarray(doses),pi=pi,R=R,LL=LL,Yk=Yk,Lk=Lk,XLk=XLk,
                prof=model.profiles(Yk,Lk,XLk))
        else:raise AssertionError('Finite decreasing-group finalization did not terminate')
        record.update(final_groups=len(current['doses']),after_R=current['R'].tolist(),after_pi=current['pi'].tolist())
        return dict(current,evidence=dict(current['evidence'],terminal_merge_checked=True,
            terminal_finalization='research_probe',search_convergence='not_established'))
    pipeline.fit_specimen=fit
    try:yield
    finally:pipeline.fit_specimen=original


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--panel-root',type=Path,required=True);p.add_argument('--assertion-source',type=Path,required=True)
    p.add_argument('--previous',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
    root=Path(__file__).parent;script=root/'joint_unordered_panel.py'
    spec=importlib.util.spec_from_file_location('locked_panel',script);panel=importlib.util.module_from_spec(spec);spec.loader.exec_module(panel)
    inputs=dict(script=Path(__file__),panel=script,previous=a.previous/'summary.json',
        previous_manifest=a.previous/'run_manifest.json',design=a.task/'DESIGN.md')
    old=json.loads(inputs['previous_manifest'].read_text())
    if old['status']!='complete':raise AssertionError('Previous screen incomplete')
    for r in old['input']['files'].values():
        if file_record(r['path'])['sha256']!=r['sha256']:raise AssertionError('Previous inputs changed')
    rec=start_run(a.out,inputs,dict(merge='margin',operation='terminal merge/refit until no merge'),
        dict(scope='Stopping-rule correctness diagnosis; no threshold or original criterion changes'),mode='terminal_merge_probe')
    trace=[];argv=sys.argv
    try:
        sys.argv=[str(script),'--out',str(a.out/'panel'),'--task',str(a.task),'--panel-root',str(a.panel_root),
            '--assertion-source',str(a.assertion_source),'--merge','margin']
        with terminal_merge(trace):panel.main()
        (a.out/'terminal_trace.json').write_text(json.dumps(trace,indent=2))
        old_cases=json.loads(inputs['previous'].read_text())['cases']
        new_cases=json.loads((a.out/'panel/summary.json').read_text())['cases']
        if [r['seed'] for r in old_cases]!=[r['seed'] for r in new_cases]:raise AssertionError('Case family changed')
        comparisons=[]
        for before,after in zip(old_cases,new_cases):
            comparisons.append(dict(seed=before['seed'],before=before['criteria'],after=after['criteria'],
                before_all_numeric=before['all_numeric_criteria_pass'],after_all_numeric=after['all_numeric_criteria_pass']))
        (a.out/'summary.json').write_text(json.dumps(comparisons,indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        (a.out/'terminal_trace_partial.json').write_text(json.dumps(trace,indent=2))
        finish_run(a.out,rec,error=error);raise
    finally:sys.argv=argv


if __name__=='__main__':main()
