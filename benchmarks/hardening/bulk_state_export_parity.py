"""Independent old-wheel capture and exact preservation audit for bulk state exports."""
import argparse
import fcntl
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

from haplodose.provenance import start_run, finish_run, file_record


CALL = r'''
import hashlib,json,sys
from pathlib import Path
import numpy as np
import haplodose, haplodose._core as native
from haplodose import bulk
data,out,mode,artifact,source,capture=sys.argv[1:]
artifact=json.loads(Path(artifact).read_text())
for key in ('native','wheel'):
    assert hashlib.sha256(Path(artifact[key]['path']).read_bytes()).hexdigest()==artifact[key]['sha256']
assert Path(native.__file__).resolve()==Path(artifact['native']['path']).resolve()
assert Path(haplodose.__file__).resolve().is_relative_to(Path(sys.executable).parent.parent)
expected=json.loads((Path(source)/'source_manifest.json').read_text())
for name,digest in expected.items():
    if name.startswith('python/haplodose/') and name.endswith(('.py','.tsv')):
        path=Path(haplodose.__file__).parent/name.removeprefix('python/haplodose/')
        assert hashlib.sha256(path.read_bytes()).hexdigest()==digest
original=bulk.fit_sample; captured=[]
def fit(*args,**kwargs):
    result=original(*args,**kwargs)
    captured.append({key:np.array(result[key],copy=True) for key in ('post','dose','path','beta')})
    return result
bulk.fit_sample=fit
data=Path(data)
bulk.run_bulk(data/'counts.tsv',data/'samples.tsv',data/'positions.tsv',out,
    bin_genes=10,min_cpm=0,n_programmes=2,bin_stat=mode,counts_kind='raw',device='cpu',
    genome_build='hg38',coordinate_system='zero-based-half-open',gene_namespace='custom')
np.savez(capture,**{key:np.stack([row[key] for row in captured]) for key in ('post','dose','path','beta')})
manifest=json.loads((Path(out)/'run_manifest.json').read_text())
assert manifest['status']=='complete' and manifest['backend']=='cpu'
print(json.dumps(dict(wheel=artifact['wheel'],fit_calls=len(captured),backend=manifest['backend'])))
'''


def contents(path):
    return gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True); parser.add_argument('--task',type=Path,required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    old_task = args.task.parent/'bin_identity_20261109'
    old_artifact = old_task/'validation_23054611/installed_artifact.json'
    new_artifact = Path(os.environ['HAPLODOSE_VALIDATED_ARTIFACT'])
    old = json.loads(old_artifact.read_text()); new = json.loads(new_artifact.read_text())
    lock = (Path(old['python']).parent.parent/'.hardening-validate.lock').open('r')
    fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
    inputs = dict(script=Path(__file__),design=args.task/'DESIGN.md',old_artifact=old_artifact,new_artifact=new_artifact,
        baseline=args.task/'baseline/source_manifest.json',candidate=args.task/'candidate_v1/source_manifest.json')
    run = start_run(args.out,inputs,dict(seed=[20261109,901],modes=['sum','median']),
        dict(scope='Public bulk preservation and independent arithmetic; no biological validation'),mode='bulk_state_export_parity')
    try:
        for artifact in (old,new):
            for key in ('native','wheel'):
                assert file_record(artifact[key]['path'])['sha256']==artifact[key]['sha256']
        baseline = json.loads(inputs['baseline'].read_text()); candidate = json.loads(inputs['candidate'].read_text())
        assert [name for name in baseline if candidate.get(name)!=baseline[name]] == ['python/haplodose/bulk.py']
        assert sorted(set(candidate)-set(baseline)) == ['tests/test_bulk_state_exports.py']
        for name,digest in baseline.items():
            assert file_record(args.task/'baseline'/name)['sha256']==digest
        data = args.out/'data'; data.mkdir()
        rng = np.random.default_rng([20261109,901]); genes = np.array([f'G{i:04d}' for i in range(120)])
        order = rng.permutation(120); names = ['T2','001','N2','T0','N0','N1','T1']
        normal = np.array([False,True,True,False,True,True,False])
        rates = np.tile(rng.lognormal(3,.8,120)[:,None],(1,7))
        rates[30:60,~normal] *= 1.35; rates[90:,~normal] *= .65
        pd.DataFrame(rng.poisson(rates*20)[order],index=genes[order],columns=names).to_csv(data/'counts.tsv',sep='\t')
        pd.DataFrame(dict(sample=names,role=np.where(normal,'normal','tumour'))).to_csv(data/'samples.tsv',sep='\t',index=False)
        within = np.tile(np.arange(60),2)
        pd.DataFrame(dict(gene=genes,chromosome=np.repeat(['chr1','chr2'],60),start=130000000+within*10000,
                          end=130000500+within*10000)).to_csv(data/'positions.tsv',sep='\t',index=False)
        (args.out/'executed_call.py').write_text(CALL)
        records = []
        for mode in ('sum','median'):
            captures = {}
            for label,artifact,source in (('old',old,args.task/'baseline'),('new',new,args.task/'candidate_v1')):
                destination = args.out/f'{mode}_{label}'; capture = args.out/f'{mode}_{label}_fits.npz'
                prefix = Path(artifact['python']).parent.parent
                env = {**os.environ,'LD_LIBRARY_PATH':str(prefix/'lib'),'JAX_PLATFORMS':'cpu'}
                result = subprocess.run([artifact['python'],'-I','-c',CALL,str(data),str(destination),mode,
                    str(old_artifact if label=='old' else new_artifact),str(source),str(capture)],
                    env=env,capture_output=True,text=True,cwd=args.out,timeout=1200)
                (args.out/f'{mode}_{label}.stdout').write_text(result.stdout)
                (args.out/f'{mode}_{label}.stderr').write_text(result.stderr)
                result.check_returncode()
                with np.load(capture) as archive: captures[label] = {key:archive[key] for key in archive.files}
                manifest = json.loads((destination/'run_manifest.json').read_text()); assert manifest['status']=='complete'
                for record in manifest['output']['files'].values():
                    assert file_record(record['path'])['sha256']==record['sha256']
            before = args.out/f'{mode}_old'; after = args.out/f'{mode}_new'
            manifest = json.loads((before/'run_manifest.json').read_text()); equal = []
            for filename in manifest['output']['files']:
                if filename == 'run.json': continue
                assert contents(before/filename)==contents(after/filename), filename
                equal.append(filename)
            previous = json.loads((before/'run.json').read_text()); current = json.loads((after/'run.json').read_text())
            metadata = current.pop('bin_rna_state_outputs'); assert previous==current
            assert metadata['calibration_status']=='unvalidated' and metadata['biological_DNA_state_probability'] is False
            for key in ('post','dose','path','beta'):
                assert np.array_equal(captures['old'][key],captures['new'][key]), key
            posterior = np.load(after/'rna_state_posterior.npy'); effect = np.load(after/'rna_log_effect.npy')
            path = np.load(after/'rna_state_index.npy')
            for key,array in (('post',posterior),('dose',effect),('path',path)):
                assert np.array_equal(array,captures['old'][key]), key
            samples = pd.read_csv(after/'samples.tsv',sep='\t',dtype={'sample':str},float_precision='round_trip')
            assert samples['sample'].tolist()==names
            beta = captures['old']['beta'][:,None,None]
            true_levels = np.log(np.arange(1,5,dtype=np.float64)/2)
            stored_levels = true_levels.astype(np.float32).astype(np.float64)
            cast_beta = beta.astype(np.float32).astype(np.float64)
            coefficients = cast_beta*stored_levels
            u = np.finfo(np.float32).eps/2; gamma = 4*u/(1-4*u)
            coefficient_error = np.abs(coefficients-beta*true_levels)+u*np.abs(coefficients)
            independent = np.sum(posterior.astype(np.float64)*beta*true_levels,axis=2)
            budget = np.sum(posterior*coefficient_error,axis=2)+gamma*np.sum(np.abs(posterior*coefficients),axis=2)
            error = np.abs(effect-independent); assert np.all(error<=budget)
            segments = pd.read_csv(after/'segments.tsv',sep='\t',dtype={'sample':str},float_precision='round_trip')
            coverage = np.zeros(effect.shape,np.int64); projected = np.zeros_like(effect)
            for row in segments.itertuples():
                i = names.index(row.sample); a,b = row.first_bin,row.last_bin
                assert np.all(path[i,a:b+1]==row.cn_state-1)
                assert row.rna_log_effect==float(effect[i,a:b+1].mean())
                coverage[i,a:b+1] += 1; projected[i,a:b+1] = row.rna_log_effect
            assert (coverage==1).all()
            change = float(np.max(np.abs(projected-effect))); assert change>0
            records.append(dict(mode=mode,historical_artifacts_exact=sorted(equal),actual_fit_arrays_exact=True,
                segments_exact=len(segments),per_bin_scalar_error_max=float(error.max()),
                per_bin_scalar_error_budget_max=float(budget.max()),maximum_segment_averaging_change=change))
        (args.out/'report.json').write_text(json.dumps(dict(modes=records,baseline_source_files_preserved=len(baseline),
            old_tests_unchanged=True,scientific_thresholds_changed=False,biological_validation=False),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error)
        raise


if __name__=='__main__': main()
