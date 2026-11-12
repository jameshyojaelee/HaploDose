"""Prespecified all-five recovery stage ablations with unchanged assertions."""
import argparse, importlib.util, json, shutil, time
from dataclasses import asdict
from pathlib import Path
import pandas as pd
from haplodose.clones import CloneParams
from haplodose.pipeline import run
from haplodose.provenance import start_run, finish_run

def load(name,path):
    sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--panel-root',type=Path,required=True);p.add_argument('--assertion-source',type=Path,required=True)
    p.add_argument('--stage',choices=['merge_control','exact_initial','exact_matched'],required=True);a=p.parse_args()
    a.out.mkdir(exist_ok=False)
    root=Path(__file__).parent
    mod=load('original_assertions',a.assertion_source)
    score=load('unchanged_scoring',root/'recovery_baseline.py')
    hook=load('weighted_stages',root/'shared_membership_stages.py')
    params=CloneParams(split_mode='thin')
    probe=load('exact_checks',root/'shared_membership_probe.py')
    model=load('exact_model',root/'shared_membership_hmm.py')
    numerical=probe.checks(model)
    import numpy as np
    depth=np.repeat(np.arange(8),6);labels=np.tile([0,0,1,1,1,1],8);ids=np.arange(48)+100
    selected,sizes=hook.balanced_subset(depth,labels,ids)
    perm=np.random.default_rng(20261108).permutation(48)
    permuted,other=hook.balanced_subset(depth[perm],labels[perm],ids[perm])
    np.testing.assert_array_equal(np.sort(ids[selected]),np.sort(ids[perm][permuted]))
    assert sizes==other and np.bincount(labels[selected]).tolist()==[16,16]
    (a.out/'checks.json').write_text(json.dumps(numerical,indent=2))
    rec=start_run(a.out,{'script':__file__,'assertions':a.assertion_source,'scorer':root/'recovery_baseline.py',
                        'experiment':root/'shared_membership_stages.py','model':root/'shared_membership_hmm.py',
                        'numerical_checks':root/'shared_membership_probe.py','merge_control':root/'membership_weighted_stages.py'},asdict(params),
                  {'stage':a.stage,'scope':'development-only, fixed criteria, no default promotion'},mode='stage_ablation')
    frozen=json.loads((a.panel_root/'manifest.json').read_text()); results=[]
    try:
        for seed in score.SEEDS:
            src=a.panel_root/f'seed_{seed}';out=a.out/f'seed_{seed}';out.mkdir();start=time.monotonic()
            entry=next(r for r in frozen['cases'] if r['seed']==seed)
            for record in entry['input_files']+[entry['model_metadata']]:
                if score.record(record['path'])['sha256']!=record['sha256']:
                    raise RuntimeError('Frozen panel input changed')
            metadata=pd.read_csv(src/'model_cells.tsv',sep='\t')
            if any(c.startswith('truth_') for c in metadata):raise RuntimeError('Truth in inference metadata')
            trace=[]
            with hook.weighted_stages(a.stage,trace):
                run(src/'data/counts.h5',src/'model_cells.tsv',src/'data/positions.tsv',out/'fit',chunk=97,params=params,**mod.INPUT_CONTRACT)
            before=pd.read_csv(src/'data/cells.tsv',sep='\t');fit=pd.read_csv(out/'fit/cells.tsv.gz',sep='\t')
            truth=[c for c in before if c.startswith('truth_')]
            joined=fit.merge(before[['cell_id']+truth],on='cell_id',validate='one_to_one');view=out/'assertion_view';view.mkdir()
            joined.to_csv(view/'cells.tsv.gz',sep='\t',index=False)
            for name in ('qc.json','differences.tsv'):shutil.copyfile(out/'fit'/name,view/name)
            ok=True
            try:mod.test_thin_records_mode_and_recovers_clones({'thin':view})
            except AssertionError:ok=False
            checks=score.criteria(joined,pd.read_csv(out/'fit/differences.tsv',sep='\t'))
            row=dict(seed=seed,original_assertion_passed=ok,criteria=checks,seconds=time.monotonic()-start,
                     s3_groups=int(joined.loc[joined.specimen.eq('S3') & joined.malignant_call,'clone'].nunique()))
            results.append(row);(out/'summary.json').write_text(json.dumps(row,indent=2))
            (out/'stage_trace.json').write_text(json.dumps(trace,indent=2));print(json.dumps(row),flush=True)
        (a.out/'summary.json').write_text(json.dumps(dict(stage=a.stage,params=asdict(params),cases=results),indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as err:
        finish_run(a.out,rec,error=err);raise

if __name__=='__main__':main()
