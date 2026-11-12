"""Strict sham then one exposed causal split-stability intervention."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def load(name):
    spec = importlib.util.spec_from_file_location('refinement_replay_'+name, Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def exact_outputs(before, after, arms=('production', 'pooled_prevalence')):
    arrays = tables = folds = 0
    for arm in arms:
        for file in sorted((before/arm/'fit').iterdir()):
            path = after/arm/'fit'/file.name
            if file.suffix == '.npy':
                np.testing.assert_array_equal(np.load(file),np.load(path)); arrays += 1
            elif file.name.endswith(('.tsv','.tsv.gz')):
                pd.testing.assert_frame_equal(pd.read_csv(file,sep='\t'),pd.read_csv(path,sep='\t'),check_exact=True); tables += 1
        for file in sorted((before/arm/'folds').glob('*.npz')):
            with np.load(file) as a, np.load(after/arm/'folds'/file.name) as b:
                assert a.files == b.files
                for key in a.files: np.testing.assert_array_equal(a[key],b[key])
            folds += 1
    return dict(arrays=arrays,tables=tables,folds=folds,exact=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--apply',action='store_true')
    args=parser.parse_args()
    trace_task=args.task.parent/'phase8_membership_trace_20261110'
    assert json.loads((trace_task/'summary_v1/run_manifest.json').read_text())['status']=='complete'
    assert json.loads((trace_task/'summary_v1/report.json').read_text())['exact_replay_gate_passed']
    if args.apply:
        assert json.loads((args.task/'sham_v1/run_manifest.json').read_text())['status']=='complete'
        assert json.loads((args.task/'sham_v1/report.json').read_text())['sham_parity']['exact']
    args.out.mkdir(exist_ok=False)
    pilot, method = load('pooled_prevalence_pilot'), load('refinement_checked_split')
    guard = pilot.load('pooled_prevalence_workflow'); original = guard.safeguard
    inputs=dict(script=Path(__file__),method=Path(method.__file__),design=args.task/'DESIGN.md',
                causal_evidence=trace_task/'summary_v1/report.json',
                original_replay=trace_task/'trace_v1/report.json')
    run=start_run(args.out,inputs,dict(seed=202611098001,apply=args.apply),
                  dict(scope='Exposed split-stability mechanism; no new error or power estimate'),
                  mode='refinement_checked_replay')
    records=[]
    @contextmanager
    def intervene(*a,**kw):
        with original(*a,**kw):
            with method.checked_splits(records,apply=args.apply): yield
    guard.safeguard=intervene; argv=sys.argv
    try:
        pilot_task=args.task.parent/'pooled_prevalence_pilot_20261110'
        sys.argv=[str(pilot.__file__),'--out',str(args.out/'replay'),'--task',str(pilot_task),
                  '--condition','two_clone','--start','1','--stop','2']
        pilot.main()
        before=trace_task/'trace_v1/replay/seed_202611098001'
        after=args.out/'replay/seed_202611098001'
        parity=exact_outputs(before,after) if not args.apply else None
        production=exact_outputs(before,after,arms=('production',)) if args.apply else None
        if args.apply:
            for file in (before/'pooled_prevalence/folds').glob('*.npz'):
                with np.load(file) as a,np.load(after/'pooled_prevalence/folds'/file.name) as b:
                    assert a.files==b.files
                    for key in a.files:np.testing.assert_array_equal(a[key],b[key])
        assert records and any(not row['passed'] for row in records)
        results=json.loads((args.out/'replay/results.json').read_text())
        write_json(dict(sham_parity=parity,production_parity=production,apply=args.apply,preview_records=records,
                        results=results,scientific_thresholds_changed=False,
                        independent_evaluation=False,biological_specimens=0,
                        promotion=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU frozen complete workflow with restored-state split previews')
    except BaseException as error:
        write_json(records,args.out/'preview_records_failed.json')
        finish_run(args.out,run,error=error); raise
    finally:
        guard.safeguard=original;sys.argv=argv


if __name__=='__main__': main()
