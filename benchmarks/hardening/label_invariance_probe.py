"""Independent full-PMF verification and three locked unordered-model diagnoses."""
import argparse
import importlib.util
import itertools
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import expit, logsumexp
from scipy.stats import nbinom
from haplodose.clones import CloneParams, _thin_diff_params
from haplodose.profile import DIFF_LEVELS, DIFF_PRIOR, credible_difference_runs
from haplodose.provenance import start_run, finish_run, file_record


def load(name):
    spec=importlib.util.spec_from_file_location('orbit_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def brute(y,E,r,p,effects,prior,change):
    """Sum full NB probabilities over labelled states and all cell assignments."""
    cells,length=y.shape;records=[];null=float(nbinom.logpmf(y,r,r/(r+E)).sum())
    for path in itertools.product(range(3),repeat=length):
        weight=prior[path[0]]
        for a,b in zip(path[:-1],path[1:]):weight*=1-change if a==b else change*prior[b]/(1-prior[a])
        oriented=list(path);nonzero=[j for j,s in enumerate(path) if s!=1]
        reverse=bool(nonzero and path[nonzero[0]]==0)
        if reverse:oriented=[2-s for s in path]
        for z in itertools.product((0,1),repeat=cells):
            means=np.array([[E[c,j]*math.exp(z[c]*effects[s])/(1-p+p*math.exp(effects[s]))
                             for j,s in enumerate(path)] for c in range(cells)])
            logweight=math.log(weight)+math.fsum(math.log(p if v else 1-p) for v in z)
            ll=float(nbinom.logpmf(y,r,r/(r+means)).sum())-null
            labels=np.array([1-v if reverse else v for v in z],float) if nonzero else np.full(cells,.5)
            records.append((tuple(oriented),labels,logweight+ll))
    evidence=float(logsumexp([x[2] for x in records]));post=np.zeros((length,3));membership=np.zeros(cells);orbits={}
    for path,labels,value in records:
        prob=math.exp(value-evidence)
        for j,s in enumerate(path):post[j,s]+=prob
        membership+=prob*labels;orbits[path]=orbits.get(path,0.)+prob
    allowance=64*np.finfo(float).eps*len(records)*max(1.,abs(null),max(abs(x[2]) for x in records))
    return dict(log_evidence=evidence,posterior=post,membership=membership,
                path=np.array(max(orbits,key=orbits.get)),allowance=allowance)


def checks(model,unordered):
    y=np.array([[2.,9.,1.],[7.,3.,6.]]);E=np.array([[3.,5.,2.],[4.,2.,5.]]);r=np.full_like(E,3.)
    effects=np.array([-.7,0.,.7]);prior=np.array([.1,.8,.1]);records=[]
    for p in (.2,.5,.8):
        m=model.build(y,E,r,np.full(len(y),p),effects,prior,.04,left_neutral=False,right_neutral=False)
        actual=unordered.enumerate_unordered(m);expected=brute(y,E,r,p,effects,prior,.04)
        for key in ('posterior','membership','log_evidence'):
            np.testing.assert_allclose(actual[key],expected[key],rtol=0,atol=expected['allowance'])
        np.testing.assert_array_equal(actual['path'],expected['path'])
        mirror=unordered.enumerate_unordered(model.build(y,E,r,np.full(len(y),1-p),effects,prior,.04,
                                             left_neutral=False,right_neutral=False))
        perm=np.array([1,0]);reordered=unordered.enumerate_unordered(model.build(y[perm],E[perm],r[perm],
            np.full(len(y),p),effects,prior,.04,left_neutral=False,right_neutral=False))
        for key in ('posterior','membership','log_evidence'):
            np.testing.assert_allclose(actual[key],mirror[key],rtol=0,atol=expected['allowance'])
            np.testing.assert_allclose(actual[key],reordered[key][perm] if key=='membership' else reordered[key],
                                       rtol=0,atol=expected['allowance'])
        labelled=model.enumerate_states(m)
        np.testing.assert_allclose(actual['posterior'][:,1],labelled['posterior'][:,1],rtol=0,atol=expected['allowance'])
        if actual['unordered_trajectories']!=(3**3+1)//2:raise AssertionError('Wrong orbit count')
        records.append(dict(p=p,independent_allowance=expected['allowance'],
            posterior_error=float(np.max(np.abs(actual['posterior']-expected['posterior']))),
            membership_error=float(np.max(np.abs(actual['membership']-expected['membership'])))))
    anchored=model.build(y,E,r,np.array([.2,.8]),effects,prior,.04,left_neutral=False,right_neutral=False)
    try:unordered.enumerate_unordered(anchored)
    except ValueError:pass
    else:raise AssertionError('Anchored memberships must not be silently relabelled')
    # Opposite signs at different loci must survive a GLOBAL label swap.
    named,flip,_,_=unordered.canonical(np.array([[0,1,2],[2,1,0],[1,1,1]]))
    np.testing.assert_array_equal(named,np.array([[2,1,0],[2,1,0],[1,1,1]]))
    return records


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--joint-probe',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
    model=load('shared_membership_hmm');unordered=load('label_invariant_membership');root=Path(__file__).parent
    old=json.loads((a.joint_probe/'run_manifest.json').read_text())
    if old['status']!='complete':raise ValueError('Joint probe incomplete')
    inputs={k:Path(r['path']) for k,r in old['input']['files'].items()}
    inputs.update(probe=root/'label_invariance_probe.py',unordered=root/'label_invariant_membership.py',
                  joint_manifest=a.joint_probe/'run_manifest.json',joint_results=a.joint_probe/'summary.json')
    for r in old['input']['files'].values():
        if file_record(r['path'])['sha256']!=r['sha256']:raise ValueError('Old input changed')
    rec=start_run(a.out,inputs,dict(estimand='unordered two-component state trajectory, globally oriented at first nonneutral bin'),
        dict(scope='Postfit exact label-invariance investigation; not independent confirmation, coverage or pipeline recovery'),mode='label_invariance')
    try:
        (a.out/'checks.json').write_text(json.dumps(checks(model,unordered),indent=2))
        previous=json.loads(inputs['joint_results'].read_text())['cases'];pp=_thin_diff_params(CloneParams(split_mode='thin').profile)
        results=[]
        for seed in (20261108,20261109,20261110):
            events=[x for x in json.loads(inputs[f'{seed}:trace.json'].read_text()) if x['stage']=='event_fit']
            if len(events)!=1:raise ValueError('Original event changed')
            with np.load(inputs[f'{seed}:event{events[0]["index"]}']) as z:
                O,T,E,r,cols,parameters=[z[k] for k in ('O','Ytest','E','r','cols','parameters')]
            chrom=pd.read_csv(inputs[f'{seed}:fit/bins.tsv'],sep='\t').chromosome.to_numpy()
            beta=float(pd.read_csv(inputs[f'{seed}:fit/specimens.tsv'],sep='\t').set_index('specimen').loc['S1','beta'])
            p=float(expit(parameters[2]));fits=[]
            for probability in (p,1-p):
                m=model.build(O+T,2*E,2*r,np.full(len(O),probability),DIFF_LEVELS*beta*np.log(1.5),DIFF_PRIOR,pp.p_change,
                    left_neutral=bool(cols[0]>0 and chrom[cols[0]-1]==chrom[cols[0]]),
                    right_neutral=bool(cols[-1]+1<len(chrom) and chrom[cols[-1]+1]==chrom[cols[-1]]))
                fits.append(unordered.enumerate_unordered(m))
            first,mirror=fits;old_fit=next(x for x in previous if x['seed']==seed and x['method']=='joint')
            eps=np.finfo(float).eps;nops=64+4*O.size;gamma=nops*eps/(1-nops*eps)
            allowance=gamma*max(1.,float(np.abs(m['ll']).sum()))
            for key in ('posterior','membership','log_evidence'):
                np.testing.assert_allclose(first[key],mirror[key],rtol=0,atol=4*allowance)
            np.testing.assert_array_equal(first['path'],mirror['path'])
            np.testing.assert_allclose(first['posterior'][:,1],np.array(old_fit['posterior'])[:,1],rtol=0,atol=4*allowance)
            if abs(first['log_evidence']-old_fit['log_evidence'])>allowance:raise AssertionError('Observable model changed')
            support=first['posterior'][np.arange(len(cols)),first['path']]
            cores=[[int(cols[x]),int(cols[y]),int(s)] for x,y,s in credible_difference_runs(first['path'],first['posterior'],
                    np.r_[True,np.zeros(len(cols)-1,bool)],pp)]
            row=dict(seed=seed,cores=cores,path=first['path'].tolist(),posterior=first['posterior'].tolist(),
                minimum_oriented_state_posterior=float(support.min()),neutral_trajectory_probability=first['neutral_trajectory_probability'],
                anchor_probability=first['anchor_probability'].tolist(),log_evidence=first['log_evidence'],
                maximum_label_swap_posterior_error=float(np.max(np.abs(first['posterior']-mirror['posterior']))),
                maximum_label_swap_membership_error=float(np.max(np.abs(first['membership']-mirror['membership']))),
                numerical_allowance=allowance,numerical_boundary_unresolved=bool(np.min(np.abs(support-pp.posterior))<=4*allowance),
                previous_ordered_cores=old_fit['cores'],label_convention=first['label_convention'],
                scientific_status='Different estimand: ordered-test failures remain; this is not pipeline recovery or biological validation')
            results.append(row);np.savez_compressed(a.out/f'membership_{seed}.npz',probability=first['membership'])
            print(json.dumps({k:v for k,v in row.items() if k not in ('posterior','anchor_probability')}),flush=True)
        (a.out/'summary.json').write_text(json.dumps(results,indent=2))
        baseline=json.loads((a.task/'baseline/source_manifest.json').read_text());core=0
        for name,digest in baseline.items():
            if file_record(a.task/'baseline'/name)['sha256']!=digest:raise AssertionError('Baseline changed')
            if name.startswith(('python/','tests/','crates/')):
                if file_record(a.task.parents[1]/name)['sha256']!=digest:raise AssertionError('Core changed')
                core+=1
        for r in rec['input']['files'].values():
            if file_record(r['path'])['sha256']!=r['sha256']:raise AssertionError('Input changed')
        (a.out/'preservation.json').write_text(json.dumps(dict(baseline_files=len(baseline),core_files=core,inputs=len(rec['input']['files'])),indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
