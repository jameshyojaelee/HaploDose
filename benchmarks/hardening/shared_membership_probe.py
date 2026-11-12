"""Exact mathematical checks and the three locked conditional confirmation probes."""
import argparse,importlib.util,itertools,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import nbinom
from haplodose.clones import CloneParams, _thin_diff_params
from haplodose.profile import DIFF_LEVELS,DIFF_PRIOR,credible_difference_runs
from haplodose.provenance import start_run,finish_run


def load():
    path=Path(__file__).with_name('shared_membership_hmm.py')
    s=importlib.util.spec_from_file_location('shared_latent',path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    return m,path


def independent_all(O,E,r,q,effects,prior,change,left,right,mixture_proportion=None):
    """Explicit joint enumeration of both cells and state paths, SciPy NB PMF."""
    n,t=O.shape;prior=np.asarray(prior,float);prior=prior/prior.sum()
    paths=list(itertools.product(range(3),repeat=t));prior_weights=[]
    def move(a,b):return 1-change if a==b else change*prior[b]/(1-prior[a])
    for path in paths:
        v=move(1,path[0]) if left else prior[path[0]]
        for a,b in zip(path[:-1],path[1:]):v*=move(a,b)
        if right:v*=move(path[-1],1)
        prior_weights.append(v)
    normalizer=sum(prior_weights);logs=[];identities=[]
    base=nbinom.logpmf(O,r,r/(r+E)).sum()
    for path,prior_weight in zip(paths,prior_weights):
        for z in itertools.product((0,1),repeat=n):
            zp=np.prod([q[c] if z[c] else 1-q[c] for c in range(n)])
            if zp==0:continue
            mu=np.empty_like(E)
            for j,state in enumerate(path):
                proportion=(sum(E[c,j]*q[c] for c in range(n))/sum(E[:,j])
                            if mixture_proportion is None else mixture_proportion[j])
                denominator=1-proportion+proportion*np.exp(effects[state])
                for c in range(n):mu[c,j]=E[c,j]*np.exp(effects[state]*z[c])/denominator
            value=nbinom.logpmf(O,r,r/(r+mu)).sum()-base+np.log(prior_weight/normalizer)+np.log(zp)
            logs.append(value);identities.append((path,z))
    evidence=float(logsumexp(logs));weights=np.exp(np.asarray(logs)-evidence)
    post=np.zeros((t,3));membership=np.zeros(n)
    for w,(path,z) in zip(weights,identities):
        for j,state in enumerate(path):post[j,state]+=w
        membership+=w*np.asarray(z)
    allowance=128*np.finfo(float).eps*max(1.,float(np.max(np.abs(logs))))
    return evidence,post,membership,allowance


def checks(m):
    O=np.array([[1.,3.,2.,4.],[3.,2.,5.,1.],[7.,3.,1.,4.]])
    E=np.array([[1.5,2.,3.,2.],[2.,2.,2.,2.],[4.,3.,2.,5.]])
    r=np.full_like(E,2.3);q=np.array([.15,.65,.8]);effects=np.array([-.3,0.,.3]);prior=np.asarray(DIFF_PRIOR,float)
    results=[]
    for left,right in itertools.product((False,True),repeat=2):
        model=m.build(O,E,r,q,effects,prior,.02,left_neutral=left,right_neutral=right)
        exact=m.enumerate_states(model);ref=independent_all(O,E,r,q,effects,prior,.02,left,right)
        for x,y in [(exact['log_evidence'],ref[0]),(exact['posterior'],ref[1]),(exact['membership'],ref[2])]:
            np.testing.assert_allclose(x,y,rtol=0.,atol=ref[3])
        vi=m.variational(model);fine=m.variational(model,tolerance=np.finfo(float).eps**(2/3))
        if vi['elbo']>exact['log_evidence']+ref[3]:raise AssertionError('ELBO exceeds exact evidence')
        if not vi['converged'] or not fine['converged']:raise AssertionError('Fixture convergence failed')
        perm=np.array([2,0,1]);permuted=m.enumerate_states(m.build(O[perm],E[perm],r[perm],q[perm],effects,prior,.02,left_neutral=left,right_neutral=right))
        swapped=m.enumerate_states(m.build(O,E,r,1-q,effects,prior,.02,left_neutral=left,right_neutral=right))
        np.testing.assert_allclose(exact['posterior'],permuted['posterior'],rtol=0.,atol=ref[3])
        np.testing.assert_allclose(exact['posterior'],swapped['posterior'][:,::-1],rtol=0.,atol=ref[3])
        results.append(dict(left=left,right=right,evidence_error=exact['log_evidence']-ref[0],roundoff_allowance=ref[3],
                            elbo_gap=exact['log_evidence']-vi['elbo'],refinement_posterior_difference=float(np.max(np.abs(vi['posterior']-fine['posterior'])))))
    pure=m.build(O,E,r,np.array([0.,1.,1.]),effects,prior,.02,left_neutral=True,right_neutral=True)
    ex=m.enumerate_states(pure);vi=m.variational(pure)
    np.testing.assert_allclose(ex['posterior'],vi['posterior'],rtol=0.,atol=128*np.finfo(float).eps)
    np.testing.assert_allclose(ex['log_evidence'],vi['elbo'],rtol=0.,atol=128*np.finfo(float).eps)
    # Removing cells cannot change either retained cell's component means.
    model=m.build(O,E,r,q,effects,prior,.02,left_neutral=True,right_neutral=True)
    ix=np.array([0,2]);pi=model['mixture_proportion']
    subset=m.build(O[ix],E[ix],r[ix],q[ix],effects,prior,.02,left_neutral=True,right_neutral=True,mixture_proportion=pi)
    np.testing.assert_array_equal(subset['ll'],model['ll'][:,ix])
    ex=m.enumerate_states(subset);ref=independent_all(O[ix],E[ix],r[ix],q[ix],effects,prior,.02,True,True,pi)
    for x,y in [(ex['log_evidence'],ref[0]),(ex['posterior'],ref[1]),(ex['membership'],ref[2])]:
        np.testing.assert_allclose(x,y,rtol=0.,atol=ref[3])
    return dict(boundary_fixtures=results,pure_membership_reduction=True,label_and_cell_permutation=True,
                subset_preserves_component_rates=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--trace-root',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);m,path=load()
    rec=start_run(a.out,{'runner':__file__,'model':path},{},{'scope':'Conditional research model; fitted q/E/r and selection uncertainty not calibrated'},mode='joint_latent_probe')
    try:
        numerical=checks(m);(a.out/'checks.json').write_text(json.dumps(numerical,indent=2))
        results=[];pp=_thin_diff_params(CloneParams(split_mode='thin').profile)
        for seed in (20261108,20261109,20261110):
            root=a.trace_root/f'trace_{seed}';bins=pd.read_csv(root/'fit/bins.tsv',sep='\t')
            specimens=pd.read_csv(root/'fit/specimens.tsv',sep='\t');beta=float(specimens.set_index('specimen').loc['S1','beta'])
            for event in json.loads((root/'trace.json').read_text()):
                if event['stage']!='event_fit':continue
                z=np.load(root/event['arrays']);cols=z['cols']
                if not np.all(np.diff(cols)==1) or bins.iloc[cols].chromosome.nunique()!=1:raise ValueError('Noncontiguous event')
                c=bins.chromosome.to_numpy();left=bool(cols[0]>0 and c[cols[0]-1]==c[cols[0]])
                right=bool(cols[-1]+1<len(c) and c[cols[-1]+1]==c[cols[-1]])
                model=m.build(z['Ytest'],z['E'],z['r'],z['q'],DIFF_LEVELS*beta*np.log(1.5),DIFF_PRIOR,pp.p_change,
                              left_neutral=left,right_neutral=right)
                exact=m.enumerate_states(model);vi=m.variational(model);fine=m.variational(model,tolerance=np.finfo(float).eps**(2/3))
                if vi['elbo']>exact['log_evidence']+vi['numerical_elbo_allowance']:raise AssertionError('Probe ELBO exceeds exact reference')
                reset=np.r_[True,np.zeros(len(cols)-1,bool)]
                def cores(fit):return [(int(cols[x]),int(cols[y]),int(s)) for x,y,s in credible_difference_runs(fit['path'],fit['posterior'],reset,pp)]
                row=dict(seed=seed,trace_index=event['index'],beta=beta,p_change=pp.p_change,min_bins=pp.min_bins,posterior_cutoff=pp.posterior,
                    left_neutral=left,right_neutral=right,states_enumerated=exact['states_enumerated'],
                    exact_log_evidence=exact['log_evidence'],variational_elbo=vi['elbo'],elbo_gap=exact['log_evidence']-vi['elbo'],
                    exact_cores=cores(exact),variational_cores=cores(vi),converged=vi['converged'],iterations=vi['iterations'],
                    posterior_approximation_difference=float(np.max(np.abs(exact['posterior']-vi['posterior']))),
                    tolerance_refinement_difference=float(np.max(np.abs(fine['posterior']-vi['posterior']))),
                    exact_posterior=exact['posterior'].tolist(),variational_posterior=vi['posterior'].tolist(),
                    start_outcomes=[{k:f[k] for k in ('elbo','converged','iterations','residual')} for f in vi['starts']],
                    interpretation='Exact conditional working-model inference, not independent confirmation or calibrated CNA evidence')
                results.append(row)
                np.savez_compressed(a.out/f'seed_{seed}_posteriors.npz',rows=z['rows'],cells=z['cells'],cols=cols,
                    exact_membership=exact['membership'],variational_membership=vi['membership'],proposal_membership=z['q'])
                print(json.dumps(row),flush=True)
        (a.out/'summary.json').write_text(json.dumps(results,indent=2));finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise

if __name__=='__main__':main()
