"""Compare RNA contribution estimates on mixed spots with fixed pure profiles."""
import argparse
from itertools import combinations
import math
from pathlib import Path

import jax
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import softmax
from scipy.stats import binom

from haplodose.clones import CloneParams, SpecimenModel, neighbour_graph
from haplodose.io import write_json
from haplodose.kernels.estep import cell_loglik
from haplodose.provenance import file_record, finish_run, start_run


SEEDS = tuple(range(2026111000, 2026111016))
DEPTHS = (750, 3000)
BASELINES = ('simplex_ls', 'neighbour_ls', 'categorical', 'spatial_categorical', 'spatial_nb')
EPS = np.finfo(float).eps


def least_squares(q, profiles):
    """Best simplex reconstruction by checking the interior and all edges."""
    d = profiles[1:] - profiles[0]
    interior = np.linalg.lstsq(d.T, (q-profiles[0]).T, rcond=None)[0].T
    w = np.column_stack([1-interior.sum(axis=1), interior])
    best = np.where((w >= 0).all(axis=1), np.sum((w@profiles-q)**2, axis=1), np.inf)
    for a, b in combinations(range(3), 2):
        direction = profiles[b]-profiles[a]
        t = np.clip((q-profiles[a])@direction/(direction@direction), 0, 1)
        line = profiles[a]+t[:,None]*direction
        loss = np.sum((line-q)**2, axis=1)
        better = loss < best
        w[better] = 0
        w[better,a], w[better,b] = 1-t[better], t[better]
        best = np.minimum(best, loss)
    return w


def mixture(q, profiles, initial):
    """Maximize the conditional multinomial likelihood, with a simplex gap."""
    if (profiles.ndim != 2 or len(profiles) != 3 or np.any(profiles <= 0)
            or not np.isfinite(profiles).all()):
        raise ValueError('Three finite positive RNA profiles required')
    base = profiles[0]
    delta = (profiles[1:] - base).T
    if np.linalg.matrix_rank(delta) != 2:
        raise ValueError('Profile contributions are not identifiable on this feature panel')
    maximum = profiles.max(axis=0)
    scale = float(np.max(np.abs(np.log(profiles/base))))
    eps_ops = (16*len(base)+64)*EPS
    roundoff = eps_ops/(1-eps_ops)
    weights = np.empty_like(initial)
    radii, gaps = np.empty(len(q)), np.empty(len(q))
    for i, y in enumerate(q):
        def objective(x):
            shift = delta@x
            return -y@np.log1p(shift/base), -(y/(base+shift))@delta

        fit = minimize(objective, initial[i,1:], jac=True, method='SLSQP',
            bounds=((0, 1), (0, 1)),
            constraints=dict(type='ineq', fun=lambda x: 1-x.sum(),
                             jac=lambda x: -np.ones(2)),
            options=dict(ftol=EPS*max(1.,scale), maxiter=200))
        # Evaluate the returned point after machine-scale feasibility correction.
        x = np.clip(fit.x,0,1)
        x[1] = min(x[1],1-x[0])
        if np.max(np.abs(x-fit.x)) > roundoff:
            raise ArithmeticError('Material simplex infeasibility')
        w = np.r_[(1-x[0])-x[1],x]
        weights[i] = w
        _, gradient = objective(x)
        if not np.isfinite(gradient).all():
            raise ArithmeticError('Finite mixture gradient required')
        gap = float(gradient@x-min(0.,*gradient))
        grad_scale = np.max((y/(w@profiles))@np.abs(delta))
        gap_upper = max(0.,gap)+roundoff*(1+grad_scale)
        lower = delta.T@((y/maximum**2)[:,None]*delta)
        curvature = np.linalg.eigvalsh(lower)[0]-roundoff*np.linalg.norm(lower,ord=2)
        gaps[i] = gap_upper
        # Convexity bounds the objective gap by the Frank-Wolfe gap. The
        # Hessian lower bound converts it into a bound on the weight error.
        radii[i] = min(2., math.sqrt(16*gap_upper/curvature)) if curvature>0 else 2.
    return weights, radii, gaps


def checks():
    profiles = (np.eye(3)+1)/4
    q = np.array([[6,5,5],[8,4,4],[0,8,8]],float)/16
    target = np.array([[.5,.25,.25],[1.,0.,0.],[0.,.5,.5]])
    w, radii, gaps = mixture(q,profiles,least_squares(q,profiles))
    assert np.all(np.abs(w-target).sum(axis=1)<=radii)
    perm = np.array([2,0,1])
    pp, qq = profiles[perm][:,perm],q[::-1][:,perm]
    reverse, other_radii, _ = mixture(qq,pp,least_squares(qq,pp))
    np.testing.assert_array_less(np.abs(reverse[::-1][:,np.argsort(perm)]-w).sum(axis=1),
                                 radii+other_radii[::-1])
    identical = profiles.copy()
    identical[2] = identical[1]
    try:
        mixture(q,identical,np.full_like(q,1/3))
    except ValueError:
        pass
    else:
        raise AssertionError('Coincident profiles should be unresolved')
    # A pure population with this transcriptional state has the exact same
    # RNA probability vector as a mixture of two distinct genomic clones.
    equivalent = (3*profiles[1]+profiles[2])/4
    np.testing.assert_array_equal(equivalent,np.array([0.,.75,.25])@profiles)
    return dict(rational_interior_vertex_and_edge=True, row_feature_profile_permutation=True,
        duplicate_profile_abstention=True, state_mixture_observation_equivalence=True,
        max_fixture_gap=float(gaps.max()))


def simulate(seed, depth, stress):
    rng = np.random.default_rng(np.random.SeedSequence([seed, depth, int(stress)]))
    base = rng.lognormal(-1,1.2,120)
    effect = np.zeros((3,len(base)))
    effect[1:,:20] = .6*np.log(1.5)
    effect[2,40:60] = .6*np.log(.5)
    profiles = base*np.exp(effect)
    profiles /= profiles.sum(axis=1,keepdims=True)
    xy = np.column_stack(np.unravel_index(np.arange(100),(10,10))).astype(float)
    normal = .1+.3*xy[:,1]/9
    right = np.clip((xy[:,0]-3)/3,0,1)
    truth = np.column_stack([normal,(1-normal)*(1-right),(1-normal)*right])
    means = truth@profiles
    sampled = means.copy()
    if stress:
        state = .6*np.log(.5)*((xy[:,1]-4.5)/4.5)
        sampled[:,40:60] *= np.exp(state[:,None])
        sampled *= rng.gamma(1/.1,.1,size=sampled.shape)
        sampled /= sampled.sum(axis=1,keepdims=True)
    counts = rng.multinomial(depth,sampled)
    return profiles, xy, truth, means, counts


def fixed_spatial(ll, graph, params):
    model = object.__new__(SpecimenModel)
    model.params, model.graph, model.R_prev = params, graph, None
    return model.posterior(np.asarray(ll,np.float32),np.full(3,1/3)), model.lam


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    params = CloneParams()
    rec = start_run(a.out,dict(driver=Path(__file__),design=a.task/'DESIGN.md'),
        dict(seeds=SEEDS,depths=DEPTHS,neighbours=params.spatial_neighbors),
        dict(scope='Conditional fixed-profile synthetic RNA contribution mapping'),
        mode='spatial_fixed_profile_mixture')
    try:
        assert jax.default_backend()=='cpu'
        write_json(checks(),a.out/'numerics.json')
        rows = []
        for seed in SEEDS:
            for depth in DEPTHS:
                for stress in (False,True):
                    profiles,xy,truth,means,counts = simulate(seed,depth,stress)
                    q = counts/depth
                    graph = neighbour_graph(xy,params.spatial_neighbors,
                        cell_ids=[f'spot_{i:03d}' for i in range(len(xy))])
                    ls = least_squares(q,profiles)
                    mixed,radii,gaps = mixture(q,profiles,ls)
                    ll = counts@np.log(profiles).T
                    spatial,ll_lambda = fixed_spatial(ll,graph,params)
                    nb = cell_loglik(counts,np.full(len(counts),depth),np.log(profiles),
                                     np.full(profiles.shape[1],.1))
                    native,nb_lambda = fixed_spatial(nb,graph,params)
                    estimates = dict(mixture=mixed,simplex_ls=ls,neighbour_ls=graph@ls,
                        categorical=softmax(ll,axis=1),spatial_categorical=spatial,spatial_nb=native)
                    assert all(np.isfinite(x).all() and np.all(x >= 0) for x in estimates.values())
                    boundary = (xy[:,0]>3)&(xy[:,0]<6)
                    for method,w in estimates.items():
                        rows.append(dict(seed=seed,depth=depth,stress=stress,method=method,
                            contribution_mae=float(np.abs(w-truth).mean()),
                            regional_profile_l1=float(np.abs(w@profiles-means).sum(axis=1).mean()),
                            boundary_mae=float(np.abs(w[boundary]-truth[boundary]).mean()),
                            mixture_weight_l1_error_bound=float(radii.mean()) if method=='mixture' else 0.,
                            max_gap=float(gaps.max()) if method=='mixture' else None,
                            lambda_multinomial=float(ll_lambda),lambda_nb=float(nb_lambda)))
                    np.savez_compressed(a.out/f'{seed}_{depth}_{int(stress)}.npz',
                        profiles=profiles,xy=xy,truth=truth,means=means,counts=counts,
                        mixture_radius=radii,**estimates)
        scores = pd.DataFrame(rows)
        scores.to_csv(a.out/'scores.tsv',sep='\t',index=False)
        contrasts = []
        family = len(DEPTHS)*2*len(BASELINES)*2
        for (depth,stress),group in scores.groupby(['depth','stress']):
            mix = group[group.method.eq('mixture')].set_index('seed')
            for baseline in BASELINES:
                old = group[group.method.eq(baseline)].set_index('seed')
                np.testing.assert_array_equal(mix.index,old.index)
                for endpoint,factor in [('contribution_mae',1/3),('regional_profile_l1',1.)]:
                    gain = old[endpoint]-mix[endpoint]
                    bound = factor*mix.mixture_weight_l1_error_bound
                    wins = int((gain>bound).sum())
                    pvalue = float(binom.sf(wins-1,len(mix),.5))
                    exact_tail = sum(math.comb(len(mix),k) for k in range(wins,len(mix)+1))/2**len(mix)
                    assert abs(pvalue-exact_tail)<=len(mix)*EPS
                    contrasts.append(dict(depth=int(depth),stress=bool(stress),baseline=baseline,
                        endpoint=endpoint,n=len(mix),mean_gain=float(gain.mean()),
                        certified_wins=wins,p_adjusted=min(1.,family*pvalue)))
        assert len(contrasts)==family
        write_json(dict(contrasts=contrasts, paired_sign_family=family,
            point_and_sign_gate=all(r['mean_gain']>0 and r['p_adjusted']<=.05 for r in contrasts),
            independent_biological_specimens=0, production_changed=False),a.out/'report.json')
        for item in rec['input']['files'].values():
            assert file_record(item['path'])['sha256']==item['sha256']
        finish_run(a.out,rec,backend='CPU mixture optimization and native spatial assignments')
    except BaseException as error:
        finish_run(a.out,rec,error=error)
        raise


if __name__=='__main__':
    main()
