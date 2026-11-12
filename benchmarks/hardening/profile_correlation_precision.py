"""Fixed-budget refinement of the frozen correlation bootstrap, no new fits."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, binom
from haplodose.provenance import start_run, finish_run, file_record


def load(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def verify(record):
    if file_record(record['path'])['sha256']!=record['sha256']:
        raise AssertionError('Frozen artifact changed')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--previous',type=Path,required=True)
    parser.add_argument('--design',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    manifest_path=args.previous/'run_manifest.json';manifest=json.loads(manifest_path.read_text())
    if manifest['status']!='complete':raise AssertionError('Original report incomplete')
    for name in ('summary.json',)+tuple(f'refit_{i:03d}.npz' for i in range(9)):
        verify(manifest['output']['files'][name])
    names=('profile.py','statistics.py','subarm','pearson_reference.py','support','split')
    for name in names:verify(manifest['input']['files'][name])
    files={name:Path(manifest['input']['files'][name]['path']) for name in names}
    inputs=dict(files,script=Path(__file__),manifest=manifest_path,summary=args.previous/'summary.json',design=args.design)
    for i in range(9):inputs[f'refit_{i}']=args.previous/f'refit_{i:03d}.npz'
    record=start_run(args.out,inputs,dict(draws=1000000,checkpoints=[20000,200000,1000000],
        scientific_family=93,mc_family=279),dict(scope='Conditional numerical precision only; no new biological units'),
        mode='profile_correlation_precision')
    try:
        base=load('fixed_profile',files['profile.py']);precision=load('fixed_pearson',files['pearson_reference.py'])
        old=json.loads(inputs['summary'].read_text());original=old['comparisons']
        base.METHODS=tuple(old['methods']);base.CONTRASTS=tuple(map(tuple,old['contrasts']))
        if len(original)!=93 or len(base.METHODS)!=17:raise AssertionError('Comparator family changed')
        rows=[r for r in original if r['endpoint'] in ('rho_window','pearson_window')]
        if len(rows)!=62:raise AssertionError('Missing correlation endpoints')
        with np.load(files['support'],allow_pickle=True) as data:
            truth=data['dna'];samples=data['samples']
        metadata=pd.read_csv(files['split'],sep='\t').set_index('sample').reindex(samples)
        ci=np.flatnonzero(metadata.split.eq('calibration'));ei=np.flatnonzero(metadata.split.eq('evaluation'))
        if len(ci)!=24 or len(ei)!=24 or metadata.patient.duplicated().any():raise AssertionError('Patient roles changed')
        cache=[]
        for i in range(9):
            with np.load(inputs[f'refit_{i}'],allow_pickle=True) as z:
                np.testing.assert_array_equal(z['samples'],samples)
                d=dict(tum=samples,D=truth,W={m:z[m] for m in base.METHODS})
            common=np.isfinite(truth)&np.logical_and.reduce([np.isfinite(d['W'][m]) for m in base.METHODS])
            scores=base.correlations(d,common)
            for m in base.METHODS:
                for j in range(len(samples)):
                    x,y=d['W'][m][common[:,j],j],truth[common[:,j],j]
                    if len(x)>5 and np.ptp(x)>0 and np.ptp(y)>0:
                        scores[m]['pearson_window'][j]=precision.reference(x,y)['value']
            cache.append(np.array([scores[r['a']][r['endpoint']]-scores[r['b']][r['endpoint']] for r in rows]))
        cache=np.asarray(cache);np.save(args.out/'paired_patient_effects.npy',cache)
        # The original common-support rule gives the same missing patients for
        # every method. Compact these shared missing values BEFORE each reduction,
        # matching the original scalar scorer instead of replacing them by zeros.
        finite=np.isfinite(cache)
        if not np.array_equal(finite,np.broadcast_to(finite[:,:1],finite.shape)):
            raise AssertionError('Shared finite-patient assumption does not hold')
        points=np.array([float(v[np.isfinite(v)].mean()) for v in cache[0][:,ei]])
        if points.shape!=(62,):raise AssertionError('Point endpoint axes changed')
        np.testing.assert_array_equal(points,[r['mean_diff'] for r in rows])
        draws=np.lib.format.open_memmap(args.out/'bootstrap_draws.npy',mode='w+',dtype=np.float64,shape=(1000000,62))
        rng=np.random.default_rng(old['params']['seed'])
        for iteration in range(1000000):
            ri=int(rng.integers(8));rng.choice(ci,len(ci),replace=True)
            chosen=rng.choice(ei,len(ei),replace=True)
            chosen=chosen[finite[ri+1,0,chosen]]
            values=np.ascontiguousarray(cache[ri+1][:,chosen])
            draws[iteration]=values.mean(axis=1) if len(chosen) else np.nan
            if iteration<20000:
                independent=np.array([v.mean() if len(v) else np.nan for v in values])
                np.testing.assert_array_equal(draws[iteration],independent)
            if iteration==19999:
                for j,row in enumerate(rows):
                    prefix=np.asarray(draws[:20000,j]);prefix=prefix[np.isfinite(prefix)]
                    delta=points[j];error=prefix-delta
                    testable=bool(np.isfinite(delta) and len(prefix)>=18000 and np.ptp(prefix)>0)
                    if testable!=row['testable'] or len(prefix)!=row['valid_bootstrap_replicates']:
                        raise AssertionError('Original prefix testability changed')
                    if testable:
                        probability=float((1+np.count_nonzero(np.abs(error)>=abs(delta)))/(len(prefix)+1))
                        if probability!=row['p']:raise AssertionError('Original prefix probability changed')
                        for name,q in (('',.025),('family_',.05/(2*93))):
                            low,high=delta-np.quantile(error,[1-q,q])
                            if low!=row[name+'ci_low'] or high!=row[name+'ci_high']:
                                raise AssertionError('Original prefix interval changed')
                print(json.dumps(dict(original_correlations_reproduced=62,draws=20000)),flush=True)
            if (iteration+1)%100000==0:print(json.dumps(dict(completed=iteration+1,total=1000000)),flush=True)
        draws.flush();prefixes=[];tail=.05/(2*279)
        for n in (20000,200000,1000000):
            updated=[];mc=[]
            for j,row in enumerate(rows):
                boot=np.asarray(draws[:n,j]);boot=boot[np.isfinite(boot)];delta=points[j]
                testable=bool(np.isfinite(delta) and len(boot)>=max(99,int(.9*n)) and np.ptp(boot)>0)
                current=dict(row,valid_bootstrap_replicates=len(boot),testable=testable,
                    p=1.,ci_low=np.nan,ci_high=np.nan,family_ci_low=np.nan,family_ci_high=np.nan)
                if testable:
                    error=boot-delta;k=int(np.count_nonzero(np.abs(error)>=abs(delta)))
                    current['p']=float((k+1)/(len(boot)+1))
                    for prefix,q in (('',.025),('family_',.05/(2*93))):
                        lo,hi=delta-np.quantile(error,[1-q,q])
                        current[prefix+'ci_low'],current[prefix+'ci_high']=float(lo),float(hi)
                    sorted_error=np.sort(error);quantiles=[]
                    for q in (.05/(2*93),1-.05/(2*93)):
                        low=int(binom.ppf(tail,len(boot),q));high=int(binom.isf(tail,len(boot),q))+1
                        quantiles.append(dict(q=q,lower_rank=low,upper_rank=high,
                            lower=None if low==0 else float(sorted_error[low-1]),
                            upper=None if high>len(boot) else float(sorted_error[high-1])))
                    mc.append(dict(endpoint=row['endpoint'],a=row['a'],b=row['b'],exceedances=k,
                        p_lower=0. if k==0 else float(beta.ppf(tail,k,len(boot)-k+1)),
                        p_upper=1. if k==len(boot) else float(beta.isf(tail,k+1,len(boot)-k)),
                        error_quantile_bounds=quantiles))
                if n==20000:
                    for key in ('testable','valid_bootstrap_replicates','p','ci_low','ci_high','family_ci_low','family_ci_high'):
                        if current[key]!=row[key] and not (isinstance(current[key],float) and np.isnan(current[key]) and np.isnan(row[key])):
                            raise AssertionError(f'Original bootstrap reproduction differs: {j}:{key}')
                updated.append(current)
            # Loss results remain at their original numerical precision; no loss
            # threshold refit or new scientific endpoint enters this calculation.
            mapping={(r['endpoint'],r['a'],r['b']):r for r in updated}
            whole=[mapping.get((r['endpoint'],r['a'],r['b']),dict(r)) for r in original]
            base.adjust_family(whole)
            prefixes.append(dict(draws=n,comparisons=whole,correlation_mc=mc,unchanged_loss_bootstrap_draws=20000))
        (args.out/'summary.json').write_text(json.dumps(dict(prefixes=prefixes,
            original_correlations_reproduced=62,scientific_family=93,mc_family=279,
            scope='Same exposed patients and eight reference draws; numerical refinement, not validated biological confidence'),indent=2))
        finish_run(args.out,record,backend='cpu')
    except Exception as error:
        finish_run(args.out,record,error=error);raise


if __name__=='__main__':main()
