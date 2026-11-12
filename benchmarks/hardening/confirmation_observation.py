"""Capture the executed binary chr12 confirmation model without altering fits."""
import argparse
from contextlib import contextmanager
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from haplodose import clones
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


@contextmanager
def capture(directory,records):
    original=clones.SpecimenModel.differences
    agglomerate_code=clones.SpecimenModel.agglomerate.__code__
    profile_code=clones.difference_segments.__code__
    def observe(model,Ya,La,XLa,Yb,Lb,XLb,**kwargs):
        parent=sys._getframe(1)
        context=parent.f_locals.get('context','')
        match=re.search(r'_event_chrom\d+_bins(\d+)-(\d+)$',context)
        if parent.f_code is not agglomerate_code or not match:
            return original(model,Ya,La,XLa,Yb,Lb,XLb,**kwargs)
        bins=pd.read_csv(Path(model.rows.Y.filename).parent/'bins.tsv',sep='\t')
        lo,hi=map(int,match.groups())
        if set(bins.iloc[lo:hi+1].chromosome)!= {'chr12'}:
            return original(model,Ya,La,XLa,Yb,Lb,XLb,**kwargs)
        assert len(Ya)==len(Yb)==1 and not kwargs
        saved=[]
        def profile(frame,event,result):
            if frame.f_code is profile_code and event=='return' and result is not None:
                loc=frame.f_locals
                saved.append({key:np.asarray(loc[key]).copy() for key in
                              ('idx','sub_reset','shift','d','v','tau','var','emit','trans','init','post','path')})
        previous=sys.getprofile();assert previous is None
        sys.setprofile(profile)
        try:result=original(model,Ya,La,XLa,Yb,Lb,XLb,**kwargs)
        finally:sys.setprofile(previous)
        assert len(saved)==1
        arrays=saved[0]
        arrays.update(Ya=Ya,La=La,XLa=XLa,Yb=Yb,Lb=Lb,XLb=XLb,
                      base=model.base,gamma=model.gamma,W=model.W,rows=model.rows.rows,reset=model.reset)
        path=directory/f'comparison_{len(records):04d}.npz'
        np.savez_compressed(path,**arrays)
        seeds=[int(part.removeprefix('seed_')) for part in Path(model.rows.Y.filename).parts if part.startswith('seed_')]
        assert len(seeds)==1
        records.append(dict(arrays=str(path),seed=seeds[0],context=context,
            source=Path(model.rows.Y.filename).as_posix(),bins=str(Path(model.rows.Y.filename).parent/'bins.tsv'),
            first_bin=lo,last_bin=hi,params=asdict(model.diff_pp),beta=float(model.beta),
            programmes=model.programmes is not None,results=result))
        return result
    clones.SpecimenModel.differences=observe
    try:yield
    finally:clones.SpecimenModel.differences=original


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    a=parser.parse_args();a.out.mkdir(exist_ok=False);(a.out/'arrays').mkdir()
    path=Path(__file__).with_name('proposal_chain.py')
    spec=importlib.util.spec_from_file_location('observed_proposal_chain',path)
    runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
    rec=start_run(a.out,dict(script=Path(__file__),design=a.task/'DESIGN.md',runner=path),{},
                  dict(scope='Executed conditional confirmation model;truth-free observation'),mode='confirmation_observation')
    argv=sys.argv;records=[]
    try:
        sys.argv=[str(path),'--out',str(a.out/'replay'),'--task',str(a.task)]
        with capture(a.out/'arrays',records):runner.main()
        write_json(dict(all_public_outputs_identical=True,comparisons=len(records),production_changed=False,
                        biological_specimens=0),a.out/'report.json')
        finish_run(a.out,rec,backend='CPU executed HMM payload capture;unchanged validated inference')
    except BaseException as error:finish_run(a.out,rec,error=error);raise
    finally:
        write_json(records,a.out/'comparisons.json');sys.argv=argv


if __name__=='__main__':main()
