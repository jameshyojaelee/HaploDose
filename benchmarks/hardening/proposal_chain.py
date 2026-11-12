"""Observe candidate generation, confirmation and merge without supplying truth."""
import argparse
from contextlib import contextmanager, nullcontext
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from haplodose import clones
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


@contextmanager
def capture(directory, records):
    names=('_fit_event_mixture','_clone_candidates','_score_candidate','agglomerate','matched_differences')
    originals={name:getattr(clones if name=='_fit_event_mixture' else clones.SpecimenModel,name) for name in names}
    models={};event_calls=[]
    event_code=clones.SpecimenModel.event_mixtures.__code__
    def model_info(model):
        if model not in models:
            seed=[int(p.removeprefix('seed_')) for p in Path(model.rows.Y.filename).parts if p.startswith('seed_')]
            assert len(seed)==1 and model.alleles is None and model.graph is None
            folder=Path(model.rows.Y.filename).parent
            bins=pd.read_csv(folder/'bins.tsv',sep='\t')
            model_id=len(models)
            models[model]=(model_id,bins)
            np.savez_compressed(directory/f'model_{model_id:02d}.npz',rows=model.rows.rows,reset=model.reset)
            records.append(dict(kind='model',model=model_id,seed=seed[0],bins=str(folder/'bins.tsv')))
        return models[model]
    def save(model,kind,arrays,**metadata):
        model_id,_=model_info(model)
        path=directory/f'record_{len(records):05d}.npz'
        np.savez_compressed(path,**arrays)
        records.append(dict(kind=kind,model=model_id,arrays=str(path),**metadata))
    def mixture(O,E,r,step,**kwargs):
        frame=sys._getframe(1)
        assert frame.f_code is event_code
        loc=frame.f_locals;model=loc['self'];model_id,bins=model_info(model)
        cols=loc['cols'][loc['part']]
        result=originals['_fit_event_mixture'](O,E,r,step,**kwargs)
        source=f'event_chrom{loc["c"]}_bins{cols[0]}-{cols[-1]}'
        arrays=dict(q=result['responsibility'],cells=loc['cells'],cols=cols)
        if set(bins.iloc[cols].chromosome)=={'chr12'}:
            arrays.update(O=O,E=E,r=r)
            if result['parameters'] is not None:
                arrays.update(parameters=result['parameters'],null_logtheta=result['null_logtheta'])
        save(model,'event_fit',arrays,source=source,call=event_calls[-1],
             gain=float(result['gain']),budget=float(result['numerical_gain_budget']))
        return result
    def candidates(model,k,cells,lr,rho):
        model_info(model);call=len([r for r in records if r['kind']=='candidate_call'])
        records.append(dict(kind='candidate_call',model=models[model][0],call=call,parent=int(k)))
        event_calls.append(call)
        try:result=originals['_clone_candidates'](model,k,cells,lr,rho)
        finally:event_calls.pop()
        for _,rows,g,source in result:
            save(model,'candidate',dict(cells=rows,g=g),source=source,call=call,parent=int(k))
        return result
    def score(model,k,cells,g,source):
        before=len(model.log)
        result=originals['_score_candidate'](model,k,cells,g,source)
        save(model,'score',dict(cells=cells,g=g),source=source,parent=int(k),
             accepted=result is not None,score=None if result is None else float(result[0]),
             members=None if result is None else result[1],log=model.log[before:])
        return result
    def agglomerate(model,*args,**kwargs):
        frame=sys._getframe(1);loc=frame.f_locals
        if frame.f_code is clones.fit_specimen.__code__:
            cells,g=loc['sub'].copy(),loc['g'].copy();context='initial'
        elif frame.f_code is originals['_score_candidate'].__code__:
            cells,g=loc['cells'].copy(),loc['g'].copy();context=kwargs.get('context','')
        else:return originals['agglomerate'](model,*args,**kwargs)
        result=originals['agglomerate'](model,*args,**kwargs)
        save(model,'agglomerate',dict(cells=cells,g=g),context=context,members=result[0])
        return result
    def matched(model,labels,pairs,*args,**kwargs):
        frame=sys._getframe(1)
        result=originals['matched_differences'](model,labels,pairs,*args,**kwargs)
        if frame.f_code is clones.SpecimenModel.merge_unsupported.__code__:
            save(model,'merge_assignment',dict(labels=labels,assignment=frame.f_locals['assignment_R'],
                 doses=frame.f_locals['doses']),pairs=pairs,
                 results=[dict(testable=r.get('testable',True),segments=r['segments']) for r in result])
        return result
    replacements=dict(_fit_event_mixture=mixture,_clone_candidates=candidates,_score_candidate=score,
                      agglomerate=agglomerate,matched_differences=matched)
    for name,value in replacements.items():setattr(clones if name=='_fit_event_mixture' else clones.SpecimenModel,name,value)
    try:yield
    finally:
        for name,value in originals.items():setattr(clones if name=='_fit_event_mixture' else clones.SpecimenModel,name,value)
        models.clear()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    a=parser.parse_args();a.out.mkdir(exist_ok=False)
    root=Path(__file__).parent
    runner=load('proposal_original_runner',root/'membership_stage_panel.py')
    hook=load('proposal_known_gene',root/'gene_thinning_hook.py')
    parity=load('proposal_public_parity',root/'split_draw_stability.py').parity
    rec=start_run(a.out,dict(script=Path(__file__),design=a.task/'DESIGN.md'),{},
                  dict(scope='Original5 cases,2 canonical laws,truth-free observational capture'),mode='proposal_chain')
    argv=sys.argv
    try:
        for arm in ('bin','gene'):
            folder=a.out/arm;folder.mkdir();(folder/'arrays').mkdir();(folder/'gene_folds').mkdir()
            trace,gene_trace=[],[]
            sys.argv=[str(runner.__file__),'--out',str(folder/'panel'),'--stage','hard',
                      '--panel-root',str(a.task.parent/'recovery_20261107/recovery_panel_v1'),
                      '--assertion-source',str(a.task.parent/'recovery_20261107/baseline/tests/test_split_mode.py')]
            try:
                with hook.gene_native(gene_trace,folder/'gene_folds') if arm=='gene' else nullcontext():
                    with capture(folder/'arrays',trace):runner.main()
            finally:write_json(trace,folder/'trace.json')
            old=a.task.parent/f'split_draw_stability_20261109/gate_v2/{arm}/default/panel'
            parity(folder/'panel',old)
        write_json(dict(all_public_outputs_identical=True,production_changed=False,biological_specimens=0),a.out/'report.json')
        finish_run(a.out,rec,backend='CPU unchanged production inference; truth-free proposal capture')
    except BaseException as error:finish_run(a.out,rec,error=error);raise
    finally:sys.argv=argv


if __name__=='__main__':main()
