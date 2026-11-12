"""Independent weighted working-NB likelihood and profile-contraction audit."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import nbinom
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def reference(Y,L,x,z,dose):
    off=x[:,None]*z['gamma'][None]+z['U']@z['W']
    means=[]
    for d in dose:
        lr=np.log(z['base'])[None]+z['offset'][None]+d[None]+off
        means.append(L[:,None]*np.exp(lr-logsumexp(lr,axis=1,keepdims=True)))
    return means


def contrast(Y,R,shape,left,right):
    value=np.longdouble(0);magnitude=np.longdouble(0)
    for k,(a,b) in enumerate(zip(left,right)):
        a,b=a.astype(np.longdouble),b.astype(np.longdouble)
        y,r=Y.astype(np.longdouble),shape.astype(np.longdouble)
        first=y*np.log(b/a)
        second=(y+r)*np.log1p((b-a)/(r+a))
        weight=R[:,k,None].astype(np.longdouble)
        value+=np.sum(weight*(first-second),dtype=np.longdouble)
        magnitude+=np.sum(weight*(abs(first)+abs(second)),dtype=np.longdouble)
    return float(value),float(magnitude)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    source=a.task/'trace_v1';observed=json.loads((source/'report.json').read_text())
    rec=start_run(a.out,dict(script=Path(__file__),design=a.task/'DESIGN.md',trace=source/'report.json'),{},
                  dict(scope='Fixed-membership unpenalized working objective only; no biological validation'),mode='production_profile_audit')
    try:
        assert observed['all_public_outputs_identical'] and observed['all30_criteria_identical']
        models={};rows=[]
        for entry in observed['records']:
            if entry['model'] not in models:
                with np.load(entry['model']) as data:
                    models[entry['model']]={k:data[k] for k in data.files}
            data=models[entry['model']];Y,L,x=data['Y'],data['L'],data['x']
            folder=source/f'replay/panel/seed_{entry["seed"]}/fit'
            bins=pd.read_csv(folder/'bins.tsv',sep='\t')
            cells=pd.read_csv(folder/'cells.tsv.gz',sep='\t').iloc[data['rows']]
            assert cells.specimen.nunique()==1
            specimen=str(cells.specimen.iloc[0])
            np.testing.assert_array_equal(Y,np.load(folder/'bin_counts.npy',mmap_mode='r')[data['rows']])
            with np.load(entry['arrays']) as z:
                R=z['R'];n,J=Y.shape
                # Dense, independent count/depth-pool reconstruction.
                budget=32*n*np.finfo(float).eps
                np.testing.assert_allclose(z['Yk'],R.T@Y,rtol=budget,atol=budget)
                np.testing.assert_allclose(z['Lk'],R.T@L,rtol=budget,atol=budget)
                moments=np.column_stack([L*x,L*x*x])
                np.testing.assert_allclose(z['XLk'],R.T@moments,rtol=budget,atol=budget*np.max(abs(moments)))
                values={name:reference(Y,L,x,z,z[name]) for name in ('before','own','after')}
                shape=1/z['phi'];pmf={}
                for name,means in values.items():
                    pmf[name]=np.column_stack([nbinom.logpmf(Y,shape,shape/(shape+mu)).sum(axis=1) for mu in means])
                    assert np.isfinite(pmf[name]).all()
                tests={}
                for name,left,right in [('own','before','own'),('shared','before','after'),('sharing','own','after')]:
                    delta,magnitude=contrast(Y,R,shape,values[left],values[right])
                    direct=float(np.sum(R*(pmf[right]-pmf[left])))
                    direct_scale=float(np.sum(R*(abs(pmf[right])+abs(pmf[left]))))
                    reference_budget=64*(J+16)*np.finfo(float).eps*max(1.,magnitude,direct_scale)
                    assert abs(delta-direct)<=reference_budget
                    # Installed kernel forms each score against component0.
                    # Include both endpoint cancellation magnitudes and the
                    # softmax/reduction lengths in its float32 accumulation guard.
                    baseline=[values[left][0]]*R.shape[1]
                    _,ml=contrast(Y,R,shape,baseline,values[left])
                    _,mr=contrast(Y,R,shape,baseline,values[right])
                    operations=4*J+8*z['W'].shape[0]+64
                    eps=np.finfo(np.float32).eps
                    gamma=operations*eps/(1-operations*eps)
                    kernel_budget=gamma*max(1.,ml+mr)+reference_budget
                    kernel=entry['kernel_delta_'+name]
                    assert abs(kernel-delta)<=kernel_budget,(entry['index'],name,kernel,delta,kernel_budget)
                    tests[name]=dict(delta=delta,scipy_delta=direct,kernel_delta=kernel,
                        reference_budget=reference_budget,kernel_budget=kernel_budget,
                        direction='increase' if delta>kernel_budget else 'decrease' if delta < -kernel_budget else 'unresolved')
                mask=bins.chromosome.eq('chr12').to_numpy() & (bins.start.to_numpy()<10000000)
                pairs=[]
                for i,j in itertools.combinations(range(1,R.shape[1]),2):
                    pairs.append(dict(a=i,b=j,**{name:float(np.max(abs(z[name][i,mask]-z[name][j,mask])))
                                                for name in ('before','own','after')}))
                rows.append(dict(seed=entry['seed'],specimen=specimen,index=entry['index'],iteration=entry['iteration'],
                    components=R.shape[1],objective=tests,chr12_prefix_pair_max_log_dose=pairs,
                    exact_all_tumour_profiles_equal=all(np.array_equal(z['after'][i],z['after'][j])
                        for i,j in itertools.combinations(range(1,R.shape[1]),2)) if R.shape[1]>2 else None))
        write_json(dict(all_checks_passed=True,profile_updates=len(rows),models=len(models),rows=rows,
            production_changed=False,scientific_criteria_changed=False,biological_specimens=0,
            interpretation='Unpenalized conditional working-NB objective, not the joint profile-prior objective; no EM-monotonicity or genomic claim'),a.out/'report.json')
        finish_run(a.out,rec,backend='CPU independent SciPy NB PMFs, long-double ratios, dense pooling and installed-kernel guards')
    except BaseException as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
