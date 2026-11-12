"""Observe six frozen whole-workflow mechanisms with exact output parity."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from haplodose import clones,pipeline
from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


CASES=(('matched',0,1),('state',0,1),('state',2,2),('state',1,4),('two_clone',0,2),('two_clone',47,3))


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


@contextmanager
def capture(trace,arrays):
    posterior,fit,agglomerate=clones.posteriors,pipeline.fit_specimen,clones.SpecimenModel.agglomerate
    def save(kind,values,**metadata):
        path=arrays/f'{len(trace):04d}_{kind}.npz';np.savez_compressed(path,**values)
        trace.append(dict(index=len(trace),stage=kind,arrays=str(path),**metadata))
    def post(*args,**kwargs):
        result=posterior(*args,**kwargs);caller=sys._getframe(1)
        if caller.f_code is clones.SpecimenModel.merge_unsupported.__code__:
            loc=caller.f_locals;m=loc['self']
            save('merge_assignment',dict(rows=m.rows.rows.copy(),full_R=loc['R'].copy(),proposal_R=result.copy(),
                doses=np.asarray(loc['doses']).copy(),pi=np.asarray(loc['pi']).copy()),
                min_cells=m.params.min_cells,log_position=len(m.log))
        return result
    def agg(self,Yg,Lg,XLg,*args,**kwargs):
        context=kwargs.get('context','initial')
        result=agglomerate(self,Yg,Lg,XLg,*args,**kwargs)
        save('agglomerate',dict(rows=self.rows.rows.copy(),Y_before=np.asarray(Yg).copy(),
            Y_after=np.asarray(result[1]).copy()),context=context,members=result[0])
        return result
    def fitted(*args,**kwargs):
        result=fit(*args,**kwargs);m=result['model']
        save('final',dict(rows=m.rows.rows.copy(),R=result['R'].copy(),doses=result['doses'].copy(),
            pi=result['pi'].copy(),LL=result['LL'].copy()),log=m.log,evidence=result['evidence'])
        return result
    clones.posteriors,pipeline.fit_specimen,clones.SpecimenModel.agglomerate=post,fitted,agg
    try:yield
    finally:clones.posteriors,pipeline.fit_specimen,clones.SpecimenModel.agglomerate=posterior,fit,agglomerate


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);root=Path(__file__).parent
    gene=load('gene_trace_hook',root/'gene_thinning_hook.py');stage=load('stage_trace_hook',root/'membership_weighted_stages.py')
    veto=load('raw_veto_hook',root/'gene_veto_weights.py')
    previous=a.task.parent/'gene_veto_fresh_20261109'
    audited=json.loads((previous/'evaluation_audit_v1/report.json').read_text());assert audited['all_checks_passed']
    original=[]
    for folder in sorted(previous.glob('chunk_*_v1')):original.extend(json.loads((folder/'results.json').read_text()))
    inputs=dict(script=Path(__file__),gene=root/'gene_thinning_hook.py',stage=root/'membership_weighted_stages.py',
        veto=root/'gene_veto_weights.py',design=a.task/'DESIGN.md',prior_audit=previous/'evaluation_audit_v1/report.json')
    selections=[]
    for condition,rep,groups in CASES:
        candidates=sorted([r for r in audited['records'] if r['condition']==condition and r['arm']=='candidate' and r['reported_groups']==groups],key=lambda r:r['rep'])
        assert candidates[0]['rep']==rep
        selected=next(r for r in original if r['condition']==condition and r['rep']==rep and r['arm']=='candidate')
        selections.append(selected)
        case=Path(selected['case'])
        for name in ('counts.h5','model_cells.tsv','positions.tsv'):inputs[f'{condition}:{rep}:{name}']=case/'input'/name
    rec=start_run(a.out,inputs,{},dict(scope='Selected failure mechanisms only; prior primary events frozen'),mode='gene_state_trace')
    summaries=[]
    try:
        for selected in selections:
            seed=selected['seed'];case=Path(selected['case']);out=a.out/selected['condition']/f'seed_{seed}';out.mkdir(parents=True)
            arrays=out/'arrays';arrays.mkdir();folds=out/'folds';folds.mkdir();vdir=out/'veto';vdir.mkdir()
            traces,gene_trace,stage_trace,veto_trace=[],[],[],[]
            with gene.gene_native(gene_trace,folds,case_seeds=(seed,)),stage.weighted_stages('both',stage_trace,record_memberships=True), \
                    veto.observe_weights(veto_trace,vdir,apply=True,weight_units='raw'),capture(traces,arrays):
                pipeline.run(case/'input/counts.h5',case/'input/model_cells.tsv',case/'input/positions.tsv',out/'fit',
                    params=clones.CloneParams(split_mode='thin'),chunk=97,counts_kind='raw',genome_build='hg38',
                    coordinate_system='one-based-closed',gene_namespace='symbol')
            old=case/'candidate/fit'
            for name in ('cells.tsv.gz','clones.tsv','differences.tsv','segments.tsv','profiles.tsv.gz'):
                pd.testing.assert_frame_equal(pd.read_csv(out/'fit'/name,sep='\t'),pd.read_csv(old/name,sep='\t'),check_exact=True)
            for name in ('posteriors.npy',):
                assert np.array_equal(np.load(out/'fit'/name),np.load(old/name),equal_nan=True)
            write_json(dict(stages=traces,gene=gene_trace,weighted=stage_trace,veto=veto_trace),out/'trace.json')
            summaries.append(dict(condition=selected['condition'],rep=selected['rep'],seed=seed,prior=str(case),out=str(out),
                all_original_arrays_and_tables_identical=True,trace=file_record(out/'trace.json')))
            write_json(summaries,a.out/'report.json')
        finish_run(a.out,rec,backend='CPU observational wrappers; original gene/raw-veto candidate')
    except BaseException as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
