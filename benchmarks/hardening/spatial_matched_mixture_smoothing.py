"""Test the fixed mixture likelihood under the same smoothing as least squares."""
import argparse
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


METHOD = 'neighbour_mixture'


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
    path = Path(__file__).with_name('spatial_fixed_profile_mixture.py')
    spec = importlib.util.spec_from_file_location('frozen_spatial_design',path)
    prior = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prior)
    parent = a.task.parent/'spatial_fixed_profile_mixture_20261110/study_v1'
    inputs = dict(script=Path(__file__),design=a.task/'DESIGN.md',parent_script=path,
        prior_scores=parent/'scores.tsv', prior_manifest=parent/'run_manifest.json')
    paths = {}
    for seed in prior.SEEDS:
        for depth in prior.DEPTHS:
            for stress in (False,True):
                key = (seed,depth,stress)
                path = parent/f'{seed}_{depth}_{int(stress)}.npz'
                paths[key] = path
                inputs[path.stem] = path
    rec = start_run(a.out,inputs,dict(method=METHOD,seeds=prior.SEEDS,depths=prior.DEPTHS),
        dict(scope='Frozen conditional spatial challenge with equal neighbour averaging'),
        mode='spatial_matched_mixture_smoothing')
    try:
        assert jax.default_backend()=='cpu'
        manifest = json.loads(inputs['prior_manifest'].read_text())
        assert manifest['status']=='complete'
        frozen = {item['path']:item['sha256'] for item in manifest['output']['files'].values()}
        for item in rec['input']['files'].values():
            if Path(item['path']).parent==parent and Path(item['path']).suffix in ('.npz','.tsv'):
                assert item['sha256']==frozen[item['path']]
        before = pd.read_csv(inputs['prior_scores'],sep='\t',float_precision='round_trip')
        methods = ('mixture',*prior.BASELINES)
        assert len(before)==len(methods)*len(paths)
        index = before.set_index(['seed','depth','stress','method'],verify_integrity=True)
        rows = []
        for (seed,depth,stress),path in paths.items():
            with np.load(path) as z:
                xy,truth,profiles,means,radius = (z[k] for k in
                    ('xy','truth','profiles','means','mixture_radius'))
                boundary = (xy[:,0]>3)&(xy[:,0]<6)
                for method in methods:
                    previous = index.loc[(seed,depth,stress,method)]
                    for endpoint,value in measure(z[method],truth,profiles,means,boundary).items():
                        assert value==previous[endpoint]
                graph = neighbour_graph(xy,CloneParams().spatial_neighbors,
                    cell_ids=[f'spot_{i:03d}' for i in range(len(xy))])
                np.testing.assert_array_equal(graph@z['simplex_ls'],z['neighbour_ls'])
                np.testing.assert_array_equal(np.asarray(graph.sum(axis=1)).ravel(),np.ones(len(xy)))
                weights = graph@z['mixture']
            contribution_bound = np.asarray(graph@radius)
            np.savez_compressed(a.out/path.name,contributions=weights,contribution_bound=contribution_bound)
            rows.append(dict(seed=seed,depth=depth,stress=stress,method=METHOD,
                **measure(weights,truth,profiles,means,boundary),
                mixture_weight_l1_error_bound=float(contribution_bound.mean())))
        added = pd.DataFrame(rows)
        all_scores = pd.concat([before,added],ignore_index=True)
        all_scores.to_csv(a.out/'scores.tsv',sep='\t',index=False)
        replay = pd.read_csv(a.out/'scores.tsv',sep='\t',float_precision='round_trip')
        pd.testing.assert_frame_equal(replay.iloc[:len(before)],before,check_exact=True)
        contrasts = []
        family = len(prior.DEPTHS)*2*len(methods)*2
        for (depth,stress),group in all_scores.groupby(['depth','stress']):
            mix = group[group.method.eq(METHOD)].set_index('seed')
            for method in methods:
                base = group[group.method.eq(method)].set_index('seed')
                np.testing.assert_array_equal(mix.index,base.index)
                for endpoint,factor in [('contribution_mae',1/3),('regional_profile_l1',1.)]:
                    gain = base[endpoint]-mix[endpoint]
                    allowance = factor*(base.mixture_weight_l1_error_bound+mix.mixture_weight_l1_error_bound)
                    wins = int((gain>allowance).sum())
                    tail = sum(math.comb(len(mix),k) for k in range(wins,len(mix)+1))/2**len(mix)
                    contrasts.append(dict(depth=int(depth),stress=bool(stress),baseline=method,
                        endpoint=endpoint,n=len(mix),mean_gain=float(gain.mean()),certified_wins=wins,
                        p_adjusted=min(1.,family*tail)))
        assert len(contrasts)==family
        write_json(dict(contrasts=contrasts, paired_sign_family=family, exact_prior_scores=len(before),
            exact_prior_endpoints=3*len(before), graph_matches_every_prior_ls=True,
            point_and_sign_gate=all(r['mean_gain']>0 and r['p_adjusted']<=.05 for r in contrasts),
            biological_specimens=0, production_changed=False),a.out/'report.json')
        for item in rec['input']['files'].values():
            assert file_record(item['path'])['sha256']==item['sha256']
        finish_run(a.out,rec,backend='CPU fixed spatial neighbour projection and paired scorer')
    except BaseException as error:
        finish_run(a.out,rec,error=error)
        raise


if __name__=='__main__':
    main()
