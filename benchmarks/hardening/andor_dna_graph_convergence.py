"""Matched Leiden iteration diagnostic; graph objective is not biological truth."""
import argparse
import importlib.util
import inspect
import json
import math
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import load_npz

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def modularity(graph,labels):
    """Weighted undirected modularity using independently summed sparse entries."""
    rows,cols=graph.nonzero();weights=np.asarray(graph[rows,cols]).ravel().astype(float)
    total=math.fsum(weights);assert total>0
    internal=math.fsum(weights[labels[rows]==labels[cols]])
    degree=np.array([math.fsum(graph.data[graph.indptr[i]:graph.indptr[i+1]]) for i in range(len(labels))])
    penalty=math.fsum((math.fsum(degree[labels==k])/total)**2 for k in np.unique(labels))
    return internal/total-penalty


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    root=Path(__file__).parent
    fixed=load('fixed_graph_reference',root/'andor_dna_fixed_graph.py')
    verified=load('verified_convergence',root/'andor_dna_stability_audit.py')
    contrast=load('convergence_pairs',root/'andor_dna_contrast_stability.py')
    inputs=dict(script=Path(__file__),fixed_helper=root/'andor_dna_fixed_graph.py',
                verification_helper=root/'andor_dna_stability_audit.py',pair_helper=root/'andor_dna_contrast_stability.py',
                design=args.task/'GRAPH_CONVERGENCE_DESIGN.md',
                fixed_run=args.task/'fixed_graph_v1/run_manifest.json',
                trace_run=args.task/'cluster_trace_v2/run_manifest.json',
                reconstruction=args.task/'reconstruction_v1/run_manifest.json')
    run=start_run(args.out,inputs,dict(iterations=[2,-1],resolution=1.,seed_repeats=20),
                  dict(scope='One fixed uncertain DNA graph; no truth/model-default selection'),mode='andor_dna_graph_convergence')
    try:
        for name in ('fixed_graph_v1','cluster_trace_v2','reconstruction_v1'):
            verified.verified_run(args.task/name)
        (args.out/'installed_leiden.py').write_text(inspect.getsource(sc.tl.leiden))
        (args.out/'installed_graph_conversion.py').write_text(inspect.getsource(sc._utils.get_igraph_from_adjacency))
        graph=load_npz(args.task/'cluster_trace_v2/corrected/connectivities.npz')
        assert (graph!=graph.T).nnz==0 and np.all(graph.diagonal()==0)
        stages=pd.read_csv(args.task/'cluster_trace_v2/corrected/stages.tsv',sep='\t',dtype={'rg':str})
        anchor=contrast.read_groups(args.task/'reconstruction_v1/corrected/SNU-601/cells.tsv')
        seeds=pd.read_csv(args.task/'fixed_graph_v1/fit_groups.tsv',sep='\t')
        ig=sc._utils.get_igraph_from_adjacency(graph,directed=False)
        # Worst-case floating accumulation budget from nnz+n operations, with
        # margin for normalized degree products and the independent igraph path.
        operations=graph.nnz+graph.shape[0];eps=np.finfo(float).eps
        bound=16*operations*eps/(1-operations*eps)
        records=[];pairs=[]
        for row in seeds.itertuples(index=False):
            rep=int(row.repeat);seed=int(row.seed)
            saved=pd.read_csv(args.task/f'fixed_graph_v1/seed_{rep:02d}.tsv',sep='\t',dtype={'rg':str})
            np.testing.assert_array_equal(saved.rg,stages.rg)
            original=fixed.partition(graph,seed)
            np.testing.assert_array_equal(original,saved.leiden.to_numpy())
            a=ad.AnnData(X=np.zeros((len(stages),1),np.float32))
            a.obsp['connectivities']=graph.copy();a.uns['neighbors']=dict(connectivities_key='connectivities')
            sc.tl.leiden(a,resolution=1.,flavor='igraph',n_iterations=-1,random_state=seed)
            converged=a.obs.leiden.astype(int).to_numpy();scores=[]
            for labels in (original,converged):
                q=modularity(graph,labels);reference=ig.modularity(labels.tolist(),weights='weight',resolution=1.,directed=False)
                assert abs(q-reference)<=bound
                scores.append(q)
            assert scores[1]>=scores[0]-bound
            sizes=pd.Series(converged).value_counts()
            groups=[f'leiden_{k}' if sizes[k]>=10 else 'unassigned' for k in converged]
            candidate=pd.Series('excluded',index=anchor.index);candidate.loc[stages.rg]=groups
            pd.DataFrame(dict(rg=stages.rg,leiden=converged,initial_group=groups)).to_csv(args.out/f'seed_{rep:02d}.tsv',sep='\t',index=False)
            for item in contrast.pair_counts(anchor.tolist(),candidate.tolist()):
                item.update(repeat=rep,seed=seed);pairs.append(item)
            records.append(dict(repeat=rep,seed=seed,two_iteration_groups=int(row.supported_groups),
                                converged_supported_groups=int((sizes>=10).sum()),two_iteration_Q=scores[0],
                                converged_Q=scores[1],Q_change=scores[1]-scores[0]))
        pd.DataFrame(records).to_csv(args.out/'matched_iterations.tsv',sep='\t',index=False)
        pd.DataFrame(pairs).to_csv(args.out/'cell_pair_counts.tsv',sep='\t',index=False)
        write_json(dict(original_twenty_seed_arrays_exact=True,numerical_reference='Independent fsum modularity vs igraph',
                        numerical_operation_budget=operations,numerical_roundoff_allowance=bound,
                        converged_groups_min=min(r['converged_supported_groups'] for r in records),
                        converged_groups_max=max(r['converged_supported_groups'] for r in records),
                        minimum_Q_change=min(r['Q_change'] for r in records),maximum_Q_change=max(r['Q_change'] for r in records),
                        biological_specimens=1,new_default=False,
                        interpretation='Matched algorithm diagnostic; increased modularity cannot validate DNA clones'),args.out/'report.json')
        write_json({str(path.relative_to(args.out)):file_record(path) for path in sorted(args.out.rglob('*'))
                    if path.is_file() and path.name not in ('run_manifest.json','output_files.json')},args.out/'output_files.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
