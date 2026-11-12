"""Matched old/new wheel runs and independent exact pool-export capture."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
from haplodose.provenance import start_run, finish_run, file_record
from haplodose.simulate import write_dataset


CALL = r'''
import hashlib,json,sys
from dataclasses import replace
from pathlib import Path
import numpy as np
import haplodose, haplodose._core as native
from haplodose import pipeline
from haplodose.clones import CloneParams,clone_order
data,out,mode,artifact,source,capture=sys.argv[1:]
artifact=json.loads(Path(artifact).read_text())
for key in ('native','wheel'):
    assert hashlib.sha256(Path(artifact[key]['path']).read_bytes()).hexdigest()==artifact[key]['sha256']
assert Path(native.__file__).resolve()==Path(artifact['native']['path']).resolve()
assert Path(haplodose.__file__).resolve().is_relative_to(Path(sys.executable).parent.parent)
expected=json.loads((Path(source)/'source_manifest.json').read_text())
for name,digest in expected.items():
    if name.startswith('python/haplodose/') and name.endswith(('.py','.tsv')):
        assert hashlib.sha256((Path(haplodose.__file__).parent/name.removeprefix('python/haplodose/')).read_bytes()).hexdigest()==digest
original=pipeline.summarise;captured={}
def observe(spec,rows,res,*args,**kwargs):
    order=clone_order(res['R'])
    for name,value in dict(R=res['R'][:,order],Y=res['Yk'][order],E=res['prof']['E'][order],
                           centre=res['prof']['centre'][order],rows=rows).items():
        captured[spec+'__'+name]=np.array(value,copy=True)
    return original(spec,rows,res,*args,**kwargs)
pipeline.summarise=observe
params=CloneParams(min_cells=1000 if mode=='empty' else 10,n_programmes=0,n_pcs=5,n_neighbors=5,
                   max_groups=5,outer_rounds=2,cycles=1,spatial_neighbors=0)
data=Path(data)
qc=pipeline.run(data/'counts.h5',data/'cells.tsv',data/'positions.tsv',out,params=params,bin_genes=5,
    counts_kind='raw',device='cpu',genome_build='hg38',coordinate_system='one-based-closed',gene_namespace='symbol')
assert qc['backend']=='cpu'
np.savez(capture,**captured)
'''


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--artifact',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    old_artifact=a.task.parent/'bulk_posterior_export_20261109/validation_23056837/installed_artifact.json'
    old=json.loads(old_artifact.read_text());new=json.loads(a.artifact.read_text())
    inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md',old_artifact=old_artifact,new_artifact=a.artifact,
                baseline=a.task/'baseline/source_manifest.json',candidate=a.task/'candidate_v1/source_manifest.json')
    run=start_run(a.out,inputs,dict(modes=['fitted','empty'],refits=4),
        dict(scope='Export/inference preservation and primary-pool audit, no accuracy claim'),mode='profile_pool_export_parity')
    try:
        baseline=json.loads(inputs['baseline'].read_text());candidate=json.loads(inputs['candidate'].read_text())
        differences=[name for name,digest in baseline.items() if candidate.get(name)!=digest]
        assert differences==['python/haplodose/pipeline.py','tests/test_pipeline.py'],differences
        data=write_dataset(a.out/'data',seed=20261109,n_chrom=4,genes_per_chrom=60,extra_genes=(),
            n_types=2,ref_cells_per_type=60,depth=1200.,
            specimens=[dict(name='A',n_cells=160,normal_frac=.1,
                clones=[dict(frac=.5,events=[('chr3',0.,1.,3)]),dict(frac=.5,events=[('chr3',0.,1.,4)])])])
        records=[]
        for mode in ('fitted','empty'):
            for label,artifact,path,source in (('old',old,old_artifact,a.task/'baseline'),('new',new,a.artifact,a.task/'candidate_v1')):
                prefix=Path(artifact['python']).parent.parent
                with (prefix/'.hardening-validate.lock').open('r') as lock:
                    fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
                    result=subprocess.run([artifact['python'],'-I','-c',CALL,str(data),str(a.out/f'{mode}_{label}'),mode,
                        str(path),str(source),str(a.out/f'{mode}_{label}_capture.npz')],cwd=a.out,
                        env={**os.environ,'LD_LIBRARY_PATH':str(prefix/'lib'),'JAX_PLATFORMS':'cpu'},
                        capture_output=True,text=True,timeout=1800)
                (a.out/f'{mode}_{label}.stdout').write_text(result.stdout)
                (a.out/f'{mode}_{label}.stderr').write_text(result.stderr)
                if result.returncode:raise RuntimeError(f'{mode}/{label} failed: {result.returncode}')
            before=a.out/f'{mode}_old';after=a.out/f'{mode}_new'
            old_profile=pd.read_csv(before/'profiles.tsv.gz',sep='\t',float_precision='round_trip')
            profile=pd.read_csv(after/'profiles.tsv.gz',sep='\t',float_precision='round_trip')
            pd.testing.assert_frame_equal(old_profile,profile[old_profile.columns],check_exact=True)
            equal=[]
            for file in before.iterdir():
                if file.name in ('run_manifest.json','qc.json','profiles.tsv.gz'):continue
                if file.name.endswith(('.tsv','.tsv.gz')):
                    pd.testing.assert_frame_equal(pd.read_csv(file,sep='\t',float_precision='round_trip'),
                        pd.read_csv(after/file.name,sep='\t',float_precision='round_trip'),check_exact=True)
                else:assert file.read_bytes()==(after/file.name).read_bytes(),file.name
                equal.append(file.name)
            with np.load(a.out/f'{mode}_old_capture.npz') as x,np.load(a.out/f'{mode}_new_capture.npz') as y:
                assert x.files==y.files
                for key in x.files:assert np.array_equal(x[key],y[key]),key
                if mode=='fitted':
                    c=pd.read_csv(after/'clones.tsv',sep='\t')
                    for spec in c.specimen.unique():
                        for k,name in enumerate(c.loc[c.specimen.eq(spec)].sort_values('column').clone):
                            p=profile.loc[profile.clone.eq(name)].sort_values('bin')
                            assert np.array_equal(p.pooled_count.to_numpy(),y[spec+'__Y'][k])
                            assert np.array_equal(p.expected_count.to_numpy(),y[spec+'__E'][k])
                            assert np.array_equal(p.log_ratio_centre.to_numpy(),np.repeat(y[spec+'__centre'][k],len(p)))
                        # Independent count conservation from the exact captured responsibilities.
                        R=y[spec+'__R'];rows=y[spec+'__rows'];Y=np.load(after/'bin_counts.npy')[rows].astype(float)
                        pooled=y[spec+'__Y'];reference=R.T@Y
                        u=np.finfo(float).eps;n=len(rows);gamma=(2*n+2)*u/(1-(2*n+2)*u)
                        assert (np.abs(pooled-reference)<=gamma*np.maximum(1,np.abs(reference))).all()
                        # A single half-count offset AFTER summing is independent
                        # of the number of components represented in the file.
                        observed=pooled.sum(axis=0);expected=y[spec+'__E'].sum(axis=0)
                        rna=np.log((observed+.5)/(expected+.5))
                        np.savez(a.out/f'{spec}_pooled_RNA_effect.npz',observed=observed,expected=expected,
                                 uncentred_log_effect=rna,centred_log_effect=rna-np.median(rna))
                else:assert len(profile)==0 and not x.files
            assert {'pooled_count','expected_count','log_ratio_centre','count_units'}<=set(profile)
            records.append(dict(mode=mode,profile_rows=len(profile),historical_fields_exact=True,
                fitted_arrays_exact=True,old_output_files_exact=sorted(equal)))
        (a.out/'report.json').write_text(json.dumps(dict(modes=records,counts_exact_from_executed_fit=True,
            production_inference_unchanged=True,scientific_criteria_unchanged=True,
            biological_accuracy_improvement='not tested'),indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out,run,error=error);raise
    finish_run(a.out,run,backend='CPU old/new isolated wheel inference')


if __name__=='__main__':main()
