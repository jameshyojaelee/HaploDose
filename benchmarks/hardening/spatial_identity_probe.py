"""Check spatial row-identity invariance without assuming clone truth."""
import argparse
from pathlib import Path

import numpy as np

from haplodose.clones import neighbour_graph
from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


def canonical_graph(xy,ids,k):
    order=np.argsort(ids,kind='stable');inverse=np.argsort(order)
    graph=neighbour_graph(xy[order],k)
    return graph[inverse][:,inverse]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    run=start_run(a.out,dict(script=Path(__file__)),dict(seed=20261109,k=[1,2,4,8]),
                  dict(scope='Computational cell-identity contract; no biological clone validation'),mode='spatial_identity_probe')
    try:
        cases=dict(duplicates=np.repeat(np.arange(20,dtype=float)[:,None]*[10.,0.],3,axis=0),
                   grid=np.column_stack(np.unravel_index(np.arange(64),(8,8))).astype(float))
        results=[]
        for name,xy in cases.items():
            n=len(xy);ids=np.array([f'cell_{i:03d}' for i in range(n)])
            permutation=np.random.default_rng(20261109).permutation(n);inverse=np.argsort(permutation)
            # Independent dense physical distances, used only for these small fixtures.
            distances=np.sqrt(((xy[:,None]-xy[None])**2).sum(axis=-1))
            np.fill_diagonal(distances,np.inf)
            for k in (1,2,4,8):
                old=neighbour_graph(xy,k)
                changed=neighbour_graph(xy[permutation],k)[inverse][:,inverse]
                candidate=canonical_graph(xy,ids,k)
                candidate_changed=canonical_graph(xy[permutation],ids[permutation],k)[inverse][:,inverse]
                assert (candidate!=candidate_changed).nnz==0
                assert np.all(candidate.diagonal()==0)
                assert np.all(np.diff(candidate.indptr)==k)
                # k-nearest geometry including any legitimate boundary ties.
                cutoff=np.sort(distances,axis=1)[:,k-1]
                for i in range(n):
                    selected=candidate.indices[candidate.indptr[i]:candidate.indptr[i+1]]
                    assert np.all(distances[i,selected]<=cutoff[i])
                results.append(dict(case=name,cells=n,k=k,old_changed_directed_entries=int((old!=changed).nnz),
                                    canonical_changed_directed_entries=int((candidate!=candidate_changed).nnz)))
        assert any(r['old_changed_directed_entries']>0 for r in results)
        write_json(dict(cases=results,old_identity_contract_fails=True,canonical_fixture_invariance=True,
                        biological_validation=False,scientific_thresholds_changed=False),a.out/'report.json')
        finish_run(a.out,run,backend='cpu')
    except Exception as error:
        finish_run(a.out,run,error=error);raise


if __name__=='__main__':main()
