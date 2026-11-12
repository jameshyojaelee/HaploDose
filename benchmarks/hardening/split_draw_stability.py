"""Count-splitting redraw diagnostic with every other inference seed unchanged."""
import argparse
from contextlib import contextmanager, nullcontext
from dataclasses import asdict
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
def split_seed(effective,directory,records):
    original_set,original_gene=clones.Rows.set_thinning,clones.thin_gene_counts
    def setup(rows,seed,size,eps=clones.THIN_EPS):
        assert seed==20261031
        original_set(rows,int(effective),size,eps)
        full=rows.take(np.arange(rows.n))
        prop=rows.take(np.arange(rows.n),'prop')
        test=rows.take(np.arange(rows.n),'test')
        assert np.array_equal(prop+test,full)
        cases=[int(p.removeprefix('seed_')) for p in Path(rows.Y.filename).parts if p.startswith('seed_')]
        assert len(cases)==1
        path=directory/f'model_{len(records):02d}.npz'
        np.savez_compressed(path,full=full,prop=prop,test=test,rows=rows.rows,size=size,
                            case=np.array(cases[0]),requested_seed=np.array(seed),effective_seed=np.array(effective),share=np.array(eps))
        records.append(dict(path=str(path),case=cases[0],requested_seed=int(seed),effective_seed=int(effective)))
    def genes(Xg,gbin,Ybin,Ytest,seed):
        assert seed==20261031
        return original_gene(Xg,gbin,Ybin,Ytest,int(effective))
    clones.Rows.set_thinning,clones.thin_gene_counts=setup,genes
    try:yield
    finally:clones.Rows.set_thinning,clones.thin_gene_counts=original_set,original_gene


def parity(left,right):
    a,b=(json.loads((p/'summary.json').read_text()) for p in (left,right))
    assert len(a['cases'])==len(b['cases'])==5
    for x,y in zip(a['cases'],b['cases']):
        assert x['seed']==y['seed'] and x['criteria']==y['criteria'] and x['original_assertion_passed']==y['original_assertion_passed']
        for name in ('cells.tsv.gz','clones.tsv','differences.tsv','segments.tsv','profiles.tsv.gz'):
            pd.testing.assert_frame_equal(pd.read_csv(left/f'seed_{x["seed"]}/fit'/name,sep='\t'),
                                          pd.read_csv(right/f'seed_{y["seed"]}/fit'/name,sep='\t'),check_exact=True)
        assert np.array_equal(np.load(left/f'seed_{x["seed"]}/fit/posteriors.npy'),
                              np.load(right/f'seed_{y["seed"]}/fit/posteriors.npy'),equal_nan=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--gate',action='store_true');p.add_argument('--arm',choices=['bin','gene'])
    p.add_argument('--gate-audit',default='gate_audit_v1/report.json')
    p.add_argument('--start',type=int,default=0);p.add_argument('--stop',type=int,default=10)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    root=Path(__file__).parent
    runner=load('redraw_original_panel',root/'membership_stage_panel.py')
    hook=load('redraw_known_gene',root/'gene_thinning_hook.py')
    bank=json.loads((a.task/'seed_bank.json').read_text());assert clones.CloneParams().seed==bank['canonical_seed']==20261031
    params=json.loads(json.dumps(asdict(clones.CloneParams(split_mode='thin'))))
    inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md',seeds=a.task/'seed_bank.json',
                runner=root/'membership_stage_panel.py',hook=root/'gene_thinning_hook.py')
    if not a.gate:inputs['gate_audit']=a.task/a.gate_audit
    rec=start_run(a.out,inputs,
        dict(gate=a.gate,arm=a.arm,start=a.start,stop=a.stop,inference_params=params),
        dict(scope='Technical count redraws conditional on5 exposed matrices, no biological replication'),mode='split_draw_stability')
    argv=sys.argv;runs=[]
    def run(arm,effective,folder,rep=None):
        folder.mkdir(parents=True)
        (folder/'folds').mkdir();(folder/'gene_folds').mkdir()
        gene_trace,fold_trace=[],[]
        sys.argv=[str(runner.__file__),'--out',str(folder/'panel'),'--stage','hard',
            '--panel-root',str(a.task.parent/'recovery_20261107/recovery_panel_v1'),
            '--assertion-source',str(a.task.parent/'recovery_20261107/baseline/tests/test_split_mode.py')]
        with hook.gene_native(gene_trace,folder/'gene_folds') if arm=='gene' else nullcontext():
            with split_seed(effective,folder/'folds',fold_trace) if effective is not None else nullcontext():
                runner.main()
        for case in json.loads((folder/'panel/summary.json').read_text())['cases']:
            qc=json.loads((folder/f'panel/seed_{case["seed"]}/fit/qc.json').read_text())
            assert all(qc['params'][key]==value for key,value in params.items())
            frozen=a.task.parent/f'exact_profile_merge_20261109/production_v1/panel/seed_{case["seed"]}/fit/qc.json'
            assert qc['params']==json.loads(frozen.read_text())['params']
        row=dict(arm=arm,rep=rep,effective_seed=effective if effective is not None else bank['canonical_seed'],
                 override_enabled=effective is not None,folder=str(folder),folds=fold_trace,gene_folds=gene_trace)
        runs.append(row);write_json(runs,a.out/'runs.json')
    try:
        if a.gate:
            for arm in ('bin','gene'):
                run(arm,None,a.out/arm/'default')
                run(arm,bank['canonical_seed'],a.out/arm/'canonical')
                parity(a.out/arm/'default/panel',a.out/arm/'canonical/panel')
            parity(a.out/'bin/default/panel',a.task.parent/'exact_profile_merge_20261109/production_v1/panel')
            write_json(dict(all_public_outputs_and_criteria_identical=True,runs=runs,production_changed=False,biological_specimens=0),a.out/'report.json')
        else:
            assert a.arm and 0<=a.start<a.stop<=10
            assert json.loads((a.task/a.gate_audit).read_text())['all_checks_passed']
            for row in bank['draws'][a.start:a.stop]:
                run(a.arm,row['seed'],a.out/f'draw_{row["rep"]:02d}',rep=row['rep'])
        finish_run(a.out,rec,backend='CPU unchanged inference/search seeds; explicit count-split seed intervention')
    except BaseException as error:finish_run(a.out,rec,error=error);raise
    finally:sys.argv=argv


if __name__=='__main__':main()
