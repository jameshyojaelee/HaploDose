"""Bounded exact-input trace of the unchanged seed108 research recovery."""
import argparse
from contextlib import ExitStack, contextmanager
from dataclasses import fields, is_dataclass
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import sparse

from haplodose import clones, pipeline, profile
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def load(name):
    spec = importlib.util.spec_from_file_location('kernel_replay_'+name, Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Observer:
    def __init__(self):
        self.arrays, self.records, self.calls = {}, [], {}
        self.entered = 0

    def encode(self, value):
        if isinstance(value, np.ndarray) or hasattr(value, '__jax_array__'):
            array = np.ascontiguousarray(value)
            if array.dtype.hasobject:
                return self.encode(array.tolist())
            digest = hashlib.sha256(str((array.shape, array.dtype.str)).encode()+array.tobytes()).hexdigest()
            if digest not in self.arrays:
                self.arrays[digest] = array.copy()
            return dict(array=digest, shape=list(array.shape), dtype=array.dtype.str)
        if sparse.issparse(value):
            matrix = value.tocsr()
            return dict(csr=self.encode((matrix.data, matrix.indices, matrix.indptr, matrix.shape)))
        if is_dataclass(value):
            return {field.name:self.encode(getattr(value, field.name)) for field in fields(value)}
        if isinstance(value, dict):
            return {str(k):self.encode(v) for k,v in value.items()}
        if isinstance(value, (list, tuple)):
            return [self.encode(v) for v in value]
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, (str, int, float, bool)) or value is None:
            return value
        raise TypeError(f'Unrecordable numerical argument: {type(value)}')

    @contextmanager
    def event_profiles(self):
        # The existing stage hook checks the fitter's DIRECT caller frame.
        # Observe its profiling callback without inserting a fitter wrapper.
        original_setprofile = sys.setprofile
        event_code = clones._fit_event_mixture.__code__
        def install(callback):
            if callback is None:
                return original_setprofile(None)
            def composed(frame, event, value):
                callback(frame, event, value)
                if frame.f_code is event_code and event == 'return' and isinstance(value, dict):
                    label = 'event_mixture'
                    self.entered += 1
                    ordinal = self.calls.get(label,0)
                    self.calls[label] = ordinal+1
                    self.records.append(dict(entered=self.entered,label=label,ordinal=ordinal,
                        inputs=self.encode({k:frame.f_locals[k] for k in ('O','E','r','step')}),
                        outputs=self.encode(value)))
            return original_setprofile(composed)
        if sys.getprofile() is not None:
            raise RuntimeError('Observer cannot replace an active profiler')
        sys.setprofile = install
        try:
            yield
        finally:
            original_setprofile(None)
            sys.setprofile = original_setprofile

    @contextmanager
    def wrap(self, owner, name, label, model=False):
        original = getattr(owner, name)
        def observed(*args, **kwargs):
            self.entered += 1
            ordinal = self.calls.get(label, 0)
            self.calls[label] = ordinal+1
            entry = dict(entered=self.entered, label=label, ordinal=ordinal)
            offset = 1 if model else 0
            entry['inputs'] = self.encode(dict(args=args[offset:], kwargs=kwargs))
            if model:
                obj = args[0]
                entry['state'] = self.encode({k:getattr(obj,k) for k in
                    ('L','base','offset','gamma','U','W','phi','reset','centre','beta','tau') if hasattr(obj,k)})
            self.records.append(entry)
            result = original(*args, **kwargs)
            entry['outputs'] = self.encode(result)
            return result
        setattr(owner, name, observed)
        try:
            yield
        finally:
            setattr(owner, name, original)


def worker(out, task):
    out.mkdir(exist_ok=False)
    case = out/'seed_20261108'
    case.mkdir()
    fit_out = case/'fit'
    panel = task.parent/'recovery_20261107/recovery_panel_v1'
    source = panel/'seed_20261108'
    frozen = json.loads((panel/'manifest.json').read_text())
    row = next(r for r in frozen['cases'] if r['seed']==20261108)
    for record in row['input_files']+[row['model_metadata']]:
        assert file_record(record['path'])['sha256'] == record['sha256']
    assertion_path = task.parent/'recovery_20261107/baseline/tests/test_split_mode.py'
    spec = importlib.util.spec_from_file_location('original_recovery_assertions', assertion_path)
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    noise, hook, stages, veto, scorer = [load(name) for name in
        ('reference_gene_noise','gene_thinning_hook','membership_weighted_stages','gene_veto_weights','recovery_baseline')]
    meta = pd.read_csv(source/'model_cells.tsv', sep='\t')
    assert not any(c.startswith('truth_') for c in meta)
    for name in ('folds','weights'):
        (out/name).mkdir()
    observer = Observer()
    estimates, splits, weights, stages_trace = [], [], [], []
    from threadpoolctl import threadpool_info
    import jax
    run = start_run(out, dict(script=Path(__file__),design=task/'DESIGN.md', assertions=assertion_path,
        noise=Path(noise.__file__),hook=Path(hook.__file__),stage=Path(stages.__file__),
        veto=Path(veto.__file__),scorer=Path(scorer.__file__),manifest=panel/'manifest.json'),
        dict(seed=20261108,chunk=97,threads=int(os.environ['OMP_NUM_THREADS'])),
        dict(scope='Unchanged single-case mechanism trace, not error/power or a passing historical exact gate'),
        mode='recovery_kernel_replay_worker')
    try:
        with ExitStack() as stack:
            for owner,name,label,model in (
                (pipeline,'reference_rates','reference_rates',False),
                (pipeline,'learn_programmes','learn_programmes',False),
                (profile,'segment','profile_segment',False),
                (clones,'segment','clones_segment',False),
                (clones.SpecimenModel,'profiles','model_profiles',True),
                (clones.SpecimenModel,'loglik','model_loglik',True)):
                stack.enter_context(observer.wrap(owner,name,label,model))
            stack.enter_context(observer.event_profiles())
            stack.enter_context(veto.observe_weights(weights,out/'weights',apply=True,weight_units='effective'))
            stack.enter_context(hook.gene_native(splits,out/'folds',gene_size=noise.estimator(estimates,5)))
            stack.enter_context(stages.weighted_stages('both',stages_trace))
            pipeline.run(source/'data/counts.h5',source/'model_cells.tsv',source/'data/positions.tsv',
                fit_out,chunk=97,params=clones.CloneParams(split_mode='thin'),**original.INPUT_CONTRACT)
        truth = pd.read_csv(source/'data/cells.tsv',sep='\t')
        fit = pd.read_csv(fit_out/'cells.tsv.gz',sep='\t')
        joined = fit.merge(truth[['cell_id']+[c for c in truth if c.startswith('truth_')]],
            on='cell_id',validate='one_to_one')
        criteria = scorer.criteria(joined,pd.read_csv(fit_out/'differences.tsv',sep='\t'))
        view = out/'assertion_view'
        view.mkdir()
        joined.to_csv(view/'cells.tsv.gz',sep='\t',index=False)
        import shutil
        for name in ('qc.json','differences.tsv'):
            shutil.copyfile(fit_out/name,view/name)
        passed = True
        try:
            original.test_thin_records_mode_and_recovers_clones({'thin':view})
        except AssertionError:
            passed = False
        np.savez_compressed(out/'arrays.npz',**observer.arrays)
        write_json(sorted(observer.records,key=lambda r:r['entered']),out/'trace.json')
        write_json(dict(criterion_records=criteria,original_assertion_passed=passed,
            calls=observer.calls,unique_arrays=len(observer.arrays),
            array_bytes=sum(a.nbytes for a in observer.arrays.values()),
            threads={k:os.environ.get(k) for k in ('OMP_NUM_THREADS','RAYON_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','PYTHONHASHSEED')},
            affinity=sorted(os.sched_getaffinity(0)),thread_pools=threadpool_info(),
            backend=jax.default_backend(),devices=list(map(str,jax.devices())),
            cpu_info=Path('/proc/cpuinfo').read_text(),
            historical_gate_passed=False),out/'report.json')
        write_json(dict(estimates=estimates,splits=splits,weights=weights,stages=stages_trace),out/'research_hooks.json')
        finish_run(out,run,backend=jax.default_backend())
    except BaseException as error:
        np.savez_compressed(out/'partial_arrays.npz',**observer.arrays)
        write_json(sorted(observer.records,key=lambda r:r['entered']),out/'partial_trace.json')
        finish_run(out,run,error=error)
        raise


def compare(left, right):
    a,b = [json.loads((p/'trace.json').read_text()) for p in (left,right)]
    by_key = lambda records:{(r['label'],r['ordinal']):r for r in records}
    aa,bb = by_key(a),by_key(b)
    diffs=[]
    with np.load(left/'arrays.npz') as ax,np.load(right/'arrays.npz') as bx:
        def first_difference(x,y):
            if x==y:
                return None
            if isinstance(x,dict) and isinstance(y,dict):
                if 'array' in x and 'array' in y:
                    xx,yy=ax[x['array']],bx[y['array']]
                    d=dict(left=x,right=y)
                    if xx.shape==yy.shape and xx.dtype.kind in 'biuf' and yy.dtype.kind in 'biuf':
                        d['unequal_elements']=int(np.sum(xx!=yy))
                        d['max_abs_difference']=float(np.max(np.abs(xx.astype(float)-yy.astype(float)))) if xx.size else 0.
                    return d
                for key in x.keys()|y.keys():
                    if x.get(key)!=y.get(key):
                        return dict(key=key,difference=first_difference(x.get(key),y.get(key)))
            if isinstance(x,list) and isinstance(y,list):
                for i,(xx,yy) in enumerate(zip(x,y)):
                    if xx!=yy:return dict(index=i,difference=first_difference(xx,yy))
            return dict(left=x,right=y)
        for key in sorted(aa.keys()&bb.keys(),key=lambda k:aa[k]['entered']):
            x,y=aa[key],bb[key]
            if x['inputs']!=y['inputs'] or x.get('state')!=y.get('state') or x['outputs']!=y['outputs']:
                diffs.append(dict(label=key[0],ordinal=key[1],entered_left=x['entered'],
                    inputs_exact=x['inputs']==y['inputs'],state_exact=x.get('state')==y.get('state'),
                    outputs_exact=x['outputs']==y['outputs'],
                    input_difference=first_difference(x['inputs'],y['inputs']),
                    output_difference=first_difference(x['outputs'],y['outputs'])))
    ar,br = [json.loads((p/'report.json').read_text()) for p in (left,right)]
    return dict(left=str(left),right=str(right),left_calls=len(a),right_calls=len(b),
        only_left=[list(k) for k in sorted(aa.keys()-bb.keys())],only_right=[list(k) for k in sorted(bb.keys()-aa.keys())],
        scientific_criteria_exact=ar['criterion_records']==br['criterion_records'],
        differing_calls=len(diffs),first_differences=diffs[:30],
        first_same_input_different_output=next((r for r in diffs if r['inputs_exact'] and r['state_exact'] and not r['outputs_exact']),None))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--worker',action='store_true')
    args=parser.parse_args()
    if args.worker:
        worker(args.out,args.task)
        return
    args.out.mkdir(exist_ok=False)
    run=start_run(args.out,dict(script=Path(__file__),design=args.task/'DESIGN.md'),
        dict(seed=20261108,thread_sequence=[2,4,2,4],pythonhashseed=0),
        dict(scope='Bounded reproducibility mechanism investigation; previous failed gate remains failed'),mode='recovery_kernel_replay')
    try:
        for i,threads in enumerate((2,4,2,4)):
            env=os.environ.copy()
            env.update(OMP_NUM_THREADS=str(threads),RAYON_NUM_THREADS=str(threads),PYTHONHASHSEED='0')
            cmd=[sys.executable,'-I',str(Path(__file__).resolve()),'--out',str(args.out/f'worker_{i}'),
                '--task',str(args.task),'--worker']
            with (args.out/f'worker_{i}.log').open('w') as stream:
                subprocess.run(cmd,env=env,stdout=stream,stderr=subprocess.STDOUT,check=True)
        comparisons=[compare(args.out/f'worker_{a}',args.out/f'worker_{b}') for a,b in ((0,2),(1,3),(0,1),(2,3))]
        write_json(dict(comparisons=comparisons,historical_gate_replaced=False,production_changed=False,
            scientific_criteria_changed=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU sequential fresh-process numerical trace')
    except BaseException as error:
        finish_run(args.out,run,error=error)
        raise


if __name__=='__main__':
    main()
