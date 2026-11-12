"""Fit mixed spots with independently estimated reference count variance."""
import argparse
from fractions import Fraction
import importlib.util
import math
from pathlib import Path

import jax
import numpy as np
import pandas as pd
from scipy.linalg import cholesky, solve_triangular
from scipy.special import softmax

from haplodose.clones import CloneParams, neighbour_graph
from haplodose.io import write_json
from haplodose.kernels.estep import cell_loglik
from haplodose.provenance import file_record, finish_run, start_run


METHOD = 'reference_gls_neighbour'
SEEDS = tuple(range(2026111020,2026111036))
LAWS = ('multinomial','state','overdispersion','state_and_overdispersion')
REFERENCE_SPOTS = 100


def load_prior():
    path = Path(__file__).with_name('spatial_fixed_profile_mixture.py')
    spec = importlib.util.spec_from_file_location('fixed_spatial_model',path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def conditional_covariance(p, depth, phi):
    """First-order covariance for normalized independent Gamma rate noise."""
    normalizer = np.eye(len(p))-p[:,None]
    sampling = (np.diag(p)-np.outer(p,p))/depth
    rate = (normalizer*p[None,:])@(normalizer*p[None,:]).T
    return sampling+(1-1/depth)*phi*rate


def reference_dispersion(reference_counts, p, depth):
    """Moment estimate using reference-only compositional variance."""
    q = reference_counts/depth
    variance = np.var(q,axis=0,ddof=1).sum()
    p2 = p@p
    sampling = (1-p2)/depth
    excess_per_phi = (1-1/depth)*np.sum(p*p*(1-2*p+p2))
    return max(0.,float((variance-sampling)/excess_per_phi))


def gls(q, profiles, depth, phi, least_squares):
    """Exact face-enumerating fit on whitened compositional coordinates."""
    p = profiles.mean(axis=0)
    covariance = conditional_covariance(p,depth,phi)
    lower = cholesky(covariance[:-1,:-1],lower=True)
    def whiten(x):
        return solve_triangular(lower,x[:,:-1].T,lower=True).T
    rates, observed = whiten(profiles),whiten(q)
    w = least_squares(observed,rates)
    delta = (rates[1:]-rates[0]).T
    hessian = delta.T@delta
    eps = np.finfo(float).eps
    budget = (16*len(p)+64)*eps
    curvature = np.linalg.eigvalsh(hessian)[0]-budget*np.linalg.norm(hessian,ord=2)
    if curvature <= 0:
        raise ArithmeticError('Unresolved whitened profile rank')
    residual = w@rates-observed
    g = residual@delta
    fw = np.sum(g*w[:,1:],axis=1)-np.minimum(0,g.min(axis=1))
    roundoff = budget*(1+np.max(np.abs(residual)@np.abs(delta),axis=1))
    radius = np.minimum(2,np.sqrt(16*(np.maximum(0,fw)+roundoff)/curvature))
    return w,radius


def checks(prior):
    frac = [Fraction(1,2),Fraction(1,4),Fraction(1,4)]
    p = np.array(frac,float)
    phi, depth = Fraction(1,8), 16
    exact = []
    for i, pi in enumerate(frac):
        row = []
        for j, pj in enumerate(frac):
            rate = sum(pk*pk*(int(i==k)-pi)*(int(j==k)-pj)
                       for k,pk in enumerate(frac))
            row.append((pi*int(i==j)-pi*pj)/depth+(1-Fraction(1,depth))*phi*rate)
        exact.append(row)
    np.testing.assert_array_equal(conditional_covariance(p,depth,float(phi)),np.asarray(exact,float))
    profiles = (np.eye(3)+1)/4
    q = np.array([[6,5,5],[8,4,4],[0,8,8]],float)/16
    expected = np.array([[.5,.25,.25],[1,0,0],[0,.5,.5]])
    w,r = gls(q,profiles,depth,float(phi),prior.least_squares)
    assert np.all(np.abs(w-expected).sum(axis=1)<=r)
    perm = np.array([1,2,0])
    reversed_w,reversed_r = gls(q[::-1,perm],profiles[:,perm],depth,float(phi),prior.least_squares)
    assert np.all(np.abs(reversed_w[::-1]-w).sum(axis=1)<=r+reversed_r[::-1])
    return dict(rational_covariance_exact=True,analytic_simplex_fixtures=True,
                omitted_bin_and_spot_permutation=True,maximum_fixture_radius=float(r.max()))


def draw(rng, depth, mean, overdispersion):
    if overdispersion:
        mean = mean*rng.gamma(10,.1,size=mean.shape)
        mean /= mean.sum(axis=1,keepdims=True)
    return rng.multinomial(depth,mean)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    prior = load_prior()
    params = CloneParams()
    rec = start_run(a.out,dict(driver=Path(__file__),design=a.task/'DESIGN.md',prior=Path(prior.__file__)),
        dict(seeds=SEEDS,depths=prior.DEPTHS,laws=LAWS,reference_spots=REFERENCE_SPOTS),
        dict(scope='Synthetic variance intervention with independent reference spots; fixed profiles'),
        mode='spatial_reference_gls')
    try:
        assert jax.default_backend()=='cpu'
        write_json(checks(prior),a.out/'numerics.json')
        rows = []
        for seed in SEEDS:
            for depth in prior.DEPTHS:
                profiles,xy,truth,means,original = prior.simulate(seed,depth,False)
                graph = neighbour_graph(xy,params.spatial_neighbors,
                    cell_ids=[f'spot_{i:03d}' for i in range(len(xy))])
                boundary = (xy[:,0]>3)&(xy[:,0]<6)
                ref = {}
                for noise in (False,True):
                    rng = np.random.default_rng(np.random.SeedSequence([seed,depth,int(noise),1]))
                    ref[noise] = draw(rng,depth,np.tile(profiles[0],(REFERENCE_SPOTS,1)),noise)
                for k,law in enumerate(LAWS):
                    state = k in (1,3)
                    noise = k in (2,3)
                    mu = means.copy()
                    if state:
                        state_effect = .6*np.log(.5)*((xy[:,1]-4.5)/4.5)
                        mu[:,40:60] *= np.exp(state_effect[:,None])
                        mu /= mu.sum(axis=1,keepdims=True)
                    rng = np.random.default_rng(np.random.SeedSequence([seed,depth,k,2]))
                    counts = original if k==0 else draw(rng,depth,mu,noise)
                    phi = reference_dispersion(ref[noise],profiles[0],depth)
                    q = counts/depth
                    ls = prior.least_squares(q,profiles)
                    mixed,radius,_ = prior.mixture(q,profiles,ls)
                    variance_mix,variance_radius = gls(q,profiles,depth,phi,prior.least_squares)
                    ll = counts@np.log(profiles).T
                    spatial,_ = prior.fixed_spatial(ll,graph,params)
                    nb = cell_loglik(counts,np.full(len(counts),depth),np.log(profiles),
                                     np.full(profiles.shape[1],.1))
                    native,_ = prior.fixed_spatial(nb,graph,params)
                    methods = dict(simplex_ls=(ls,0),neighbour_ls=(graph@ls,0),
                        categorical=(softmax(ll,axis=1),0),spatial_categorical=(spatial,0),spatial_nb=(native,0),
                        mixture=(mixed,radius.mean()),neighbour_mixture=(graph@mixed,(graph@radius).mean()))
                    methods[METHOD] = (graph@variance_mix,(graph@variance_radius).mean())
                    for name,(w,allowance) in methods.items():
                        assert np.isfinite(w).all() and (w>=0).all()
                        rows.append(dict(seed=seed,depth=depth,law=law,method=name,reference_phi=phi,
                            contribution_mae=float(np.abs(w-truth).mean()),
                            regional_profile_l1=float(np.abs(w@profiles-means).sum(axis=1).mean()),
                            boundary_mae=float(np.abs(w[boundary]-truth[boundary]).mean()),
                            numeric_allowance=float(allowance)))
                    np.savez_compressed(a.out/f'{seed}_{depth}_{law}.npz',profiles=profiles,xy=xy,
                        true_contributions=truth,latent_profile=means,counts=counts,reference_counts=ref[noise],
                        **{name:w for name,(w,_) in methods.items()})
        scores = pd.DataFrame(rows)
        scores.to_csv(a.out/'scores.tsv',sep='\t',index=False)
        metrics = ['contribution_mae','regional_profile_l1','boundary_mae','reference_phi']
        scores.groupby(['depth','law','method'])[metrics].mean().to_csv(a.out/'means.tsv',sep='\t')
        contrasts = []
        baselines = (*prior.BASELINES,'mixture','neighbour_mixture')
        family = len(prior.DEPTHS)*len(LAWS)*len(baselines)*2
        for (depth,law),group in scores.groupby(['depth','law']):
            new = group[group.method.eq(METHOD)].set_index('seed')
            for name in baselines:
                old = group[group.method.eq(name)].set_index('seed')
                np.testing.assert_array_equal(new.index,old.index)
                for endpoint,factor in [('contribution_mae',1/3),('regional_profile_l1',1.)]:
                    gain = old[endpoint]-new[endpoint]
                    bounds = factor*(old.numeric_allowance+new.numeric_allowance)
                    wins = int((gain>bounds).sum())
                    tail = sum(math.comb(len(new),j) for j in range(wins,len(new)+1))/2**len(new)
                    contrasts.append(dict(depth=int(depth),law=law,baseline=name,endpoint=endpoint,
                        mean_gain=float(gain.mean()),n=len(new),certified_wins=wins,p_adjusted=min(1.,family*tail)))
        assert len(contrasts)==family
        write_json(dict(contrasts=contrasts,paired_sign_family=family,simulated_specimens=len(rows)//8,
            point_and_sign_gate=all(r['mean_gain']>0 and r['p_adjusted']<=.05 for r in contrasts),
            biological_specimens=0,production_changed=False),a.out/'report.json')
        for item in rec['input']['files'].values():
            assert file_record(item['path'])['sha256']==item['sha256']
        finish_run(a.out,rec,backend='CPU independent reference dispersion and spatial composition estimation')
    except BaseException as error:
        finish_run(a.out,rec,error=error)
        raise


if __name__=='__main__':
    main()
