"""Test a spatial mean correction from independent, registered reference RNA."""
import argparse
from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path

import jax
import numpy as np
import pandas as pd

from haplodose.clones import CloneParams, neighbour_graph
from haplodose.io import write_json
from haplodose.provenance import file_record, finish_run, start_run


METHOD = 'registered_reference_gls'
MISMATCH = 'reversed_reference_gls'
BLOCK_SIZE = 20


def load_variance():
    path = Path(__file__).with_name('spatial_reference_gls.py')
    spec = importlib.util.spec_from_file_location('registered_variance_model',path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def regional_multiplier(counts, base, block_size):
    """Relative composition change from independent reference counts."""
    if counts.ndim!=2 or counts.shape[1]!=len(base) or len(base)%block_size:
        raise ValueError('Complete reference count and genomic block correspondence required')
    observed = counts.reshape(len(counts),-1,block_size).sum(axis=2)
    expected = base.reshape(-1,block_size).sum(axis=1)
    if (observed<=0).any() or (expected<=0).any():
        raise ValueError('Every registered band and block requires reference support')
    proportion = observed/observed.sum(axis=1,keepdims=True)
    return np.repeat(proportion/expected,block_size,axis=1)


def component_transform(profiles, multiplier):
    shifted = profiles*multiplier
    mass = shifted.sum(axis=1)
    return shifted/mass[:,None],mass


def latent_contributions(observed, mass):
    latent = observed/mass
    return latent/latent.sum(axis=1,keepdims=True)


def checks():
    reference = np.array([.25,.25,.125,.125,.125,.125])
    counts = np.array([[3,3,2,2,3,3]])
    expected = np.array([Fraction(3,4)]*2+[Fraction(1)]*2+[Fraction(3,2)]*2)
    np.testing.assert_array_equal(regional_multiplier(counts,reference,2)[0],expected.astype(float))
    # Direct rational mixing establishes the observed weights under a common multiplier.
    profiles = np.array([[.5,.25,.25],[.25,.5,.25],[.25,.25,.5]])
    multiplier = np.array([.5,1.,2.])
    p,mass = component_transform(profiles,multiplier)
    latent = np.array([[.25,.5,.25]])
    observed = latent*mass/(latent@mass)[:,None]
    expected_observed = np.array([[Fraction(8,37),Fraction(18,37),Fraction(11,37)]],float)
    # Each term below is a dyadic unnormalized rate; division is the sole rounding.
    direct = (latent@profiles)*multiplier
    direct /= direct.sum(axis=1,keepdims=True)
    budget = (16*profiles.size+64)*np.finfo(float).eps
    np.testing.assert_allclose(observed,expected_observed,rtol=0,atol=budget)
    np.testing.assert_allclose(observed@p,direct,rtol=0,atol=budget)
    np.testing.assert_allclose(latent_contributions(observed,mass),latent,rtol=0,atol=budget)
    assert not np.array_equal(observed,latent)
    identity,one = component_transform(profiles,np.ones(profiles.shape[1]))
    np.testing.assert_array_equal(identity,profiles)
    np.testing.assert_array_equal(latent_contributions(latent,one),latent)
    return dict(rational_block_multipliers=True,profile_mixture_transform=True,
                contribution_units_round_trip=True,identity_control_exact=True,roundoff_budget=budget)


def reference_counts(prior, profiles, xy, depth, seed, law):
    """Generate an independent reference on the query's registered coordinates."""
    pure = np.tile(profiles[0],(len(xy),1))
    if law in ('state','state_and_overdispersion'):
        effect = .6*np.log(.5)*((xy[:,1]-4.5)/4.5)
        pure[:,40:60] *= np.exp(effect[:,None])
        pure /= pure.sum(axis=1,keepdims=True)
    rng = np.random.default_rng(np.random.SeedSequence([seed,depth,prior.LAWS.index(law),3]))
    return prior.draw(rng,depth,pure,law in ('overdispersion','state_and_overdispersion'))


def fit(q, profiles, bands, multiplier, depth, phi, prior, ls):
    weights = np.empty((len(q),len(profiles)))
    allowance = np.empty(len(q))
    for band,m in enumerate(multiplier):
        selected = bands==band
        shifted,mass = component_transform(profiles,m)
        observed,radius = prior.gls(q[selected],shifted,depth,phi,ls)
        weights[selected] = latent_contributions(observed,mass)
        # L1 perturbation through positive component-wise scaling and renormalization.
        allowance[selected] = np.minimum(2,2*mass.max()/mass.min()*radius)
    return weights,allowance


def measure(w, truth, profile, means, boundary):
    return dict(contribution_mae=float(np.abs(w-truth).mean()),
                regional_profile_l1=float(np.abs(w@profile-means).sum(axis=1).mean()),
                boundary_mae=float(np.abs(w[boundary]-truth[boundary]).mean()))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    variance = load_variance()
    spatial = variance.load_prior()
    parent = a.task.parent/'spatial_reference_gls_20261110/study_v1'
    cases = {(seed,depth,law):parent/f'{seed}_{depth}_{law}.npz'
             for seed in variance.SEEDS for depth in spatial.DEPTHS for law in variance.LAWS}
    inputs = dict(driver=Path(__file__),design=a.task/'DESIGN.md',variance=Path(variance.__file__),
        scorer=Path(spatial.__file__),manifest=parent/'run_manifest.json',old_scores=parent/'scores.tsv',
        **{path.stem:path for path in cases.values()})
    rec = start_run(a.out,inputs,dict(method=METHOD,negative_control=MISMATCH,block_size=BLOCK_SIZE),
        dict(scope='Frozen synthetic query cases with independently drawn spatial reference counts'),
        mode='spatial_registered_reference')
    try:
        assert jax.default_backend()=='cpu'
        write_json(checks(),a.out/'numerics.json')
        manifest = json.loads(inputs['manifest'].read_text())
        assert manifest['status']=='complete'
        old_digests = {v['path']:v['sha256'] for v in manifest['output']['files'].values()}
        for key in ('old_scores',*(path.stem for path in cases.values())):
            item = rec['input']['files'][key]
            assert item['sha256']==old_digests[item['path']]
        calibration = {}
        # This stage has access to profile and coordinates only, never query counts or truth.
        for (seed,depth,law),path in cases.items():
            with np.load(path) as z:
                profiles,xy = z['profiles'],z['xy']
            bands = xy[:,1].astype(int)
            np.testing.assert_array_equal(bands,xy[:,1])
            counts = reference_counts(variance,profiles,xy,depth,seed,law)
            pooled = np.array([counts[bands==b].sum(axis=0) for b in range(bands.max()+1)])
            multiplier = regional_multiplier(pooled,profiles[0],BLOCK_SIZE)
            target = a.out/f'reference_{path.name}'
            np.savez_compressed(target,counts=counts,pooled=pooled,multiplier=multiplier,bands=bands)
            calibration[path.stem] = file_record(target)
        write_json(calibration,a.out/'REFERENCES_LOCKED.json')
        before = pd.read_csv(inputs['old_scores'],sep='\t',float_precision='round_trip')
        index = before.set_index(['seed','depth','law','method'],verify_integrity=True)
        baselines = (*spatial.BASELINES,'mixture','neighbour_mixture',variance.METHOD)
        assert len(before)==len(baselines)*len(cases)
        rows = []
        for (seed,depth,law),path in cases.items():
            locked = calibration[path.stem]
            assert file_record(locked['path'])==locked
            with np.load(locked['path']) as z:
                multiplier,bands = z['multiplier'],z['bands']
            with np.load(path) as z:
                xy,profiles,truth,means,counts,reference = (z[k] for k in
                    ('xy','profiles','true_contributions','latent_profile','counts','reference_counts'))
                boundary = (xy[:,0]>3)&(xy[:,0]<6)
                for method in baselines:
                    previous = index.loc[(seed,depth,law,method)]
                    for key,value in measure(z[method],truth,profiles,means,boundary).items():
                        assert value==previous[key]
            phi = variance.reference_dispersion(reference,profiles[0],depth)
            assert phi==index.loc[(seed,depth,law,variance.METHOD),'reference_phi']
            graph = neighbour_graph(xy,CloneParams().spatial_neighbors,
                cell_ids=[f'spot_{i:03d}' for i in range(len(xy))])
            predictions = {}
            for name,m in ((METHOD,multiplier),(MISMATCH,multiplier[::-1])):
                w,allowance = fit(counts/depth,profiles,bands,m,depth,phi,variance,spatial.least_squares)
                w = graph@w
                predictions[name] = w
                rows.append(dict(seed=seed,depth=depth,law=law,method=name,reference_phi=phi,
                    **measure(w,truth,profiles,means,boundary),numeric_allowance=float((graph@allowance).mean())))
            np.savez_compressed(a.out/f'predictions_{path.name}',**predictions)
        scores = pd.concat([before,pd.DataFrame(rows)],ignore_index=True)
        scores.to_csv(a.out/'scores.tsv',sep='\t',index=False)
        replay = pd.read_csv(a.out/'scores.tsv',sep='\t',float_precision='round_trip')
        pd.testing.assert_frame_equal(replay.iloc[:len(before)],before,check_exact=True)
        metrics = ['contribution_mae','regional_profile_l1','boundary_mae']
        scores.groupby(['depth','law','method'])[metrics].mean().to_csv(a.out/'means.tsv',sep='\t')
        contrasts = []
        family = len(spatial.DEPTHS)*len(variance.LAWS)*len(baselines)*2
        for (depth,law),group in scores.groupby(['depth','law']):
            candidate = group[group.method.eq(METHOD)].set_index('seed')
            for name in baselines:
                old = group[group.method.eq(name)].set_index('seed')
                np.testing.assert_array_equal(candidate.index,old.index)
                for endpoint,factor in [('contribution_mae',1/3),('regional_profile_l1',1.)]:
                    gain = old[endpoint]-candidate[endpoint]
                    bound = factor*(old.numeric_allowance+candidate.numeric_allowance)
                    wins = int((gain>bound).sum())
                    n = len(candidate)
                    tail = sum(math.comb(n,k) for k in range(wins,n+1))/2**n
                    contrasts.append(dict(depth=int(depth),law=law,baseline=name,endpoint=endpoint,
                        n=n,mean_gain=float(gain.mean()),certified_wins=wins,p_adjusted=min(1.,family*tail)))
        assert len(contrasts)==family
        write_json(dict(method=METHOD,negative_control=MISMATCH,contrasts=contrasts,
            paired_sign_family=family,original_score_rows_exact=len(before),exact_original_endpoints=3*len(before),
            reference_lock_cases=len(calibration),point_and_sign_gate=all(
                r['mean_gain']>0 and r['p_adjusted']<=.05 for r in contrasts),
            biological_specimens=0,production_changed=False),a.out/'report.json')
        for item in rec['input']['files'].values():
            assert file_record(item['path'])['sha256']==item['sha256']
        finish_run(a.out,rec,backend='CPU externally registered reference and spatial simplex GLS')
    except BaseException as error:
        finish_run(a.out,rec,error=error)
        raise


if __name__=='__main__':
    main()
