"""Replay one exposed phase8 failure with immutable method and exact output gates."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from haplodose import clones
from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


@contextmanager
def observe(records,directory):
    assert sys.getprofile() is None
    codes={clones.SpecimenModel.em.__code__:'em',
           clones.SpecimenModel.merge_unsupported.__code__:'merge'}
    pending={}
    def profile(frame,event,result):
        kind=codes.get(frame.f_code)
        if kind is None:return
        if event=='call':
            loc=frame.f_locals;model=loc['self']
            arrays=dict(rows=model.rows.rows.copy(),doses_before=np.asarray(loc['doses']).copy(),
                pi_before=np.asarray(loc['pi']).copy())
            if kind=='merge':arrays['R_before']=loc['R'].copy()
            for name in('base','offset','phi','gamma','U'):
                arrays[name+'_before']=np.asarray(getattr(model,name)).copy()
            pending[id(frame)]=(arrays,len(model.log))
        elif event=='return':
            arrays,position=pending.pop(id(frame));model=frame.f_locals['self']
            if result is None:raise RuntimeError('Instrumented inference stage failed')
            doses,pi=result[:2];arrays.update(doses_after=np.asarray(doses).copy(),pi_after=pi.copy())
            metadata=dict(stage=kind,log_start=position,log_end=len(model.log),
                log=model.log[position:],groups_before=len(arrays['doses_before'])-1,
                groups_after=len(doses)-1)
            if kind=='em':
                arrays.update(R_after=result[2].copy(),LL_after=result[3].copy(),
                    pooled_after=result[4].copy(),exposure_after=result[5].copy())
            else:
                arrays['R_after']=frame.f_locals['R'].copy()
                metadata.update(merged=bool(result[2]),
                    statistical_merge=bool(frame.f_locals['statistical_merge']))
            for name in('base','offset','phi','gamma','U'):
                arrays[name+'_after']=np.asarray(getattr(model,name)).copy()
            path=directory/f'{len(records):04d}_{kind}.npz'
            np.savez_compressed(path,**arrays)
            records.append(dict(arrays=str(path),**metadata))
            write_json(records,directory/'records.json')
    sys.setprofile(profile)
    try:yield
    finally:sys.setprofile(None)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--replicate',type=int,choices=(1,2,4),default=1)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    seed=202611098000+args.replicate
    pilot_task=args.task.parent/'pooled_prevalence_pilot_20261110'
    source=Path(__file__).with_name('pooled_prevalence_pilot.py')
    spec=importlib.util.spec_from_file_location('trace_frozen_phase8',source)
    pilot=importlib.util.module_from_spec(spec);spec.loader.exec_module(pilot)
    guard=pilot.load('pooled_prevalence_workflow');original=guard.safeguard
    directory=args.out/'stages';directory.mkdir();records=[]
    rec=start_run(args.out,dict(script=Path(__file__),design=args.task/'DESIGN.md',
        pilot=source,original_summary=pilot_task/'summary_v1/report.json'),
        dict(phase=8,replicate=args.replicate,seed=seed,arms=['production','pooled_prevalence']),
        dict(scope='Exposed mechanism replay only; no new error estimate or correction'),
        mode='phase8_membership_trace')
    @contextmanager
    def instrument(*a,**kw):
        # Install only after the weighted-stage context has entered. Preserve
        # method code identities used by its caller checks; no method wrapper.
        with original(*a,**kw):
            with observe(records,directory):yield
    guard.safeguard=instrument;argv=sys.argv
    try:
        sys.argv=[str(source),'--out',str(args.out/'replay'),'--task',str(pilot_task),
            '--condition','two_clone','--start',str(args.replicate),'--stop',str(args.replicate+1)]
        pilot.main()
        old_case=pilot_task/'pilot_two_clone_00_v1'/f'seed_{seed}'
        new_case=args.out/'replay'/f'seed_{seed}'
        parity=[]
        for arm in ('production','pooled_prevalence'):
            old,new=old_case/arm,new_case/arm;arrays=tables=folds=0
            for file in sorted((old/'fit').iterdir()):
                if file.suffix=='.npy':
                    np.testing.assert_array_equal(np.load(file),np.load(new/'fit'/file.name));arrays+=1
                elif file.name.endswith(('.tsv','.tsv.gz')):
                    pd.testing.assert_frame_equal(pd.read_csv(file,sep='\t'),
                        pd.read_csv(new/'fit'/file.name,sep='\t'),check_exact=True);tables+=1
            for file in sorted((old/'folds').glob('*.npz')):
                with np.load(file)as a,np.load(new/'folds'/file.name)as b:
                    assert a.files==b.files
                    for key in a.files:np.testing.assert_array_equal(a[key],b[key])
                folds+=1
            parity.append(dict(arm=arm,arrays=arrays,tables=tables,folds=folds,exact=True))
        assert records and any(r['stage']=='em'for r in records)
        assert any(r['stage']=='merge'for r in records)
        if args.replicate in (1,2):
            assert any(r['stage']=='merge'and r['merged']for r in records)
        write_json(dict(exact_output_parity=parity,stages=records,
            method_changed=False,scientific_criteria_changed=False,
            new_error_rate_estimate=False,biological_specimens=0),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU frozen complete phase8 replay with observational stage capture')
    except BaseException as error:finish_run(args.out,rec,error=error);raise
    finally:guard.safeguard=original;sys.argv=argv


if __name__=='__main__':main()
