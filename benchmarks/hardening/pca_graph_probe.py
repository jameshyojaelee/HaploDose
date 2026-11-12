"""Separate PCA precision, neighbour membership and graph insertion order."""
import argparse
import json
from pathlib import Path
import random

import igraph as ig
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from haplodose import clones
from haplodose.provenance import start_run,finish_run,file_record


def undirected(edges):
    edges=np.sort(np.asarray(edges),axis=1)
    return np.unique(edges[edges[:,0]!=edges[:,1]],axis=0)


def group(X,edges,params):
    graph=ig.Graph(n=len(X),edges=edges,directed=False).simplify()
    ig.set_random_number_generator(random.Random(params.seed))
    part=graph.community_leiden(objective_function='modularity',resolution=params.leiden_resolution,n_iterations=2)
    return clones.consolidate(np.asarray(part.membership),X,params),np.asarray(graph.get_edgelist())


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True);p.add_argument('--host',choices=['intel','amd'],required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    inputs={'runner':__file__};saved={}
    for replicate in (150,167):
        saved[replicate]={}
        for host in ('intel','amd'):
            path=a.task/f'{host}_v2/replicate_{replicate}/summary.json'
            entries=json.loads(path.read_text())['trace']
            first=next(r for r in entries if r['stage']=='leiden_groups')
            inputs[f'{replicate}:{host}:summary']=path
            inputs[f'{replicate}:{host}:arrays']=Path(first['arrays']['path'])
            if file_record(first['arrays']['path'])['sha256']!=first['arrays']['sha256']:raise ValueError('Frozen graph changed')
            with np.load(first['arrays']['path']) as z:saved[replicate][host]={name:z[name].copy() for name in z.files}
    rec=start_run(a.out,inputs,dict(host=a.host),dict(scope='Two frozen initial graphs; no inference change'),mode='pca_graph_probe')
    try:
        params=clones.CloneParams(n_programmes=0,spatial_neighbors=0,seed=20261107);reports=[]
        for replicate,data in saved.items():
            np.testing.assert_array_equal(data['intel']['F'],data['amd']['F'])
            F=data['intel']['F'];n=len(F);k=min(params.n_neighbors,n-1)
            old_graphs=[set(map(tuple,undirected(data[h]['edges']))) for h in ('intel','amd')]
            report=dict(replicate=replicate,old_undirected_edges=[len(g) for g in old_graphs],
                        old_graph_symmetric_difference=len(old_graphs[0]^old_graphs[1]),variants=[])
            for precision in ('float32','float64'):
                values=F.astype(precision)
                X=PCA(n_components=min(params.n_pcs,n-1,F.shape[1]-1),svd_solver='randomized',
                      random_state=params.seed).fit_transform(values)
                dist,idx=NearestNeighbors(n_neighbors=k+1).fit(X).kneighbors(X)
                edges=np.column_stack([np.repeat(np.arange(n),k),idx[:,1:].ravel()])
                canonical=undirected(edges)
                g,graph=group(X,edges,params);ordered,ordered_graph=group(X,canonical,params)
                if precision=='float32':
                    for name,actual in [('X',X),('idx',idx),('edges',edges),('groups',g)]:
                        np.testing.assert_array_equal(actual,data[a.host][name])
                # Independent double distance calculation. Self is removed by
                # identity; exact ties use row identity, not sorting accidents.
                squared=cdist(X.astype(float),X.astype(float),metric='sqeuclidean')
                np.fill_diagonal(squared,np.inf)
                reference=np.array([np.lexsort((np.arange(n),row))[:k] for row in squared])
                edge_ref=undirected(np.column_stack([np.repeat(np.arange(n),k),reference.ravel()]))
                current_set=set(map(tuple,canonical));reference_set=set(map(tuple,edge_ref))
                ordered_distance=np.sort(squared,axis=1)
                gap=ordered_distance[:,k]-ordered_distance[:,k-1]
                # Standard gamma_d scale for accumulation of nonnegative
                # squared distances; it excludes PCA approximation uncertainty.
                eps=np.finfo(X.dtype).eps;gamma=X.shape[1]*eps/(1-X.shape[1]*eps)
                allowance=gamma*(ordered_distance[:,k]+ordered_distance[:,k-1])
                path=a.out/f'{replicate}_{precision}.npz'
                np.savez_compressed(path,X=X,idx=idx,edges=canonical,groups=g,canonical_groups=ordered,
                                    distance_reference_edges=edge_ref,cutoff_gap=gap,distance_rounding_scale=allowance)
                report['variants'].append(dict(precision=precision,arrays=file_record(path),
                    graph_edge_count=len(canonical),edge_difference_vs_distance_reference=len(current_set^reference_set),
                    self_not_first=int(np.count_nonzero(idx[:,0]!=np.arange(n))),
                    self_edges_in_legacy_selection=int(np.count_nonzero(idx[:,1:]==np.arange(n)[:,None])),
                    graph_insertion_order_changes_membership=bool(not np.array_equal(g,ordered)),
                    gaps_within_distance_rounding_scale=int(np.count_nonzero(gap<=allowance)),
                    minimum_cutoff_squared_distance_gap=float(gap.min())))
            reports.append(report)
        # Minimal identity counterexample: equal features need not return self
        # first. This fixture diagnoses the legacy skip-first assumption only.
        duplicate=np.zeros((6,2));_,idx=NearestNeighbors(n_neighbors=3).fit(duplicate).kneighbors(duplicate)
        after=idx[:,1:];safe=np.array([row[row!=i][:2] for i,row in enumerate(idx)])
        if any(i in row for i,row in enumerate(safe)):raise AssertionError('Explicit self exclusion failed')
        (a.out/'summary.json').write_text(json.dumps(dict(cases=reports,duplicate_feature_fixture=dict(
            legacy_self_edges=int(np.count_nonzero(after==np.arange(6)[:,None])),
            explicit_identity_exclusion_self_edges=0),scope='Graph mechanism only; no biological validation or default promotion'),indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
