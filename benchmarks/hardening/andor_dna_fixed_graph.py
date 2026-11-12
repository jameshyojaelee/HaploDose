"""Separate Leiden seed sensitivity from sequencing changes on a saved DNA graph."""
import argparse
import importlib.util
import json
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


def partition(graph,seed):
    a=ad.AnnData(X=np.zeros((graph.shape[0],1),np.float32))
    a.obsp['connectivities']=graph.copy()
    a.uns['neighbors']=dict(connectivities_key='connectivities')
    sc.tl.leiden(a,resolution=1.,flavor='igraph',n_iterations=2,random_state=seed)
    return a.obs.leiden.astype(int).to_numpy()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    root=Path(__file__).parent
    verified=load('verified_graph',root/'andor_dna_stability_audit.py')
    contrast=load('graph_pairs',root/'andor_dna_contrast_stability.py')
    inputs=dict(script=Path(__file__),verified_helper=root/'andor_dna_stability_audit.py',
                pair_helper=root/'andor_dna_contrast_stability.py',
                design=args.task/'FIXED_GRAPH_DESIGN.md')
    for name in ('cluster_trace_v2','half_stages_v1','contrasts_v1','reconstruction_v1'):
        inputs[name]=args.task/name/'run_manifest.json'
    run=start_run(args.out,inputs,dict(seed_entropy=[20261109,1402],repeats=20),
                  dict(scope='Exploratory algorithmic stability of one fixed uncertain DNA graph'),mode='andor_dna_fixed_graph')
    try:
        for name in ('cluster_trace_v2','half_stages_v1','contrasts_v1'):
            verified.verified_run(args.task/name)
        reconstruction=verified.verified_run(args.task/'reconstruction_v1')
        record=reconstruction['input']['files']['DNA_source'];source=Path(record['path'])
        assert file_record(source)==record
        module=load('executed_graph_DNA',source)
        directory=args.task/'cluster_trace_v2/corrected'
        graph=load_npz(directory/'connectivities.npz')
        stages=pd.read_csv(directory/'stages.tsv',sep='\t',dtype={'rg':str},keep_default_na=False)
        assert stages.rg.is_unique and graph.shape==(len(stages),len(stages))
        graph_before=graph.copy()
        np.testing.assert_array_equal(partition(graph,module.SEED),stages.initial_leiden.to_numpy())
        # Independently reconcile all4800 final-stage pair records with the prior
        # collector, and verify the claimed absence of merge-induced changes.
        half=args.task/'half_stages_v1'
        a=pd.read_csv(half/'stage_pair_counts.tsv',sep='\t')
        b=pd.read_csv(args.task/'contrasts_v1/cell_pair_counts.tsv',sep='\t')
        keys=['pair','half','anchor_left','anchor_right']
        columns=['fixed_reference_pairs','both_assigned_pairs','coassigned_pairs','separated_pairs','unresolved_pairs']
        before=a.loc[a.stage=='initial'].set_index(keys).sort_index()
        after=a.loc[a.stage=='final'].set_index(keys).sort_index()
        prior=b.set_index(keys).sort_index()
        assert before.index.is_unique and len(before)==4800
        pd.testing.assert_index_equal(after.index,prior.index)
        np.testing.assert_array_equal(before[columns],after[columns])
        np.testing.assert_array_equal(after[columns],prior[columns])
        anchor=contrast.read_groups(args.task/'reconstruction_v1/corrected/SNU-601/cells.tsv')
        assert set(stages.rg)<=set(anchor.index)
        rows=[];fits=[];seeds=[]
        for rep in range(20):
            seed=int(np.random.SeedSequence([20261109,1402,rep]).generate_state(1)[0]);seeds.append(seed)
            labels=partition(graph,seed);sizes=pd.Series(labels).value_counts()
            groups=[f'leiden_{k}' if sizes[k]>=10 else 'unassigned' for k in labels]
            candidate=pd.Series('excluded',index=anchor.index)
            candidate.loc[stages.rg]=groups
            pd.DataFrame(dict(rg=stages.rg,leiden=labels,initial_group=groups)).to_csv(args.out/f'seed_{rep:02d}.tsv',sep='\t',index=False)
            records=contrast.pair_counts(anchor.tolist(),candidate.tolist())
            for item in records:
                item.update(repeat=rep,seed=seed)
            rows.extend(records)
            fits.append(dict(repeat=rep,seed=seed,initial_groups=len(sizes),
                             supported_groups=int((sizes>=10).sum()),unassigned_cells=int(candidate.eq('unassigned').sum())))
        assert len(set(seeds))==20
        assert (graph!=graph_before).nnz==0
        table=pd.DataFrame(rows);table.to_csv(args.out/'cell_pair_counts.tsv',sep='\t',index=False)
        pd.DataFrame(fits).to_csv(args.out/'fit_groups.tsv',sep='\t',index=False)
        summaries=[]
        for (left,right),group in table.groupby(['anchor_left','anchor_right'],sort=True):
            assert len(group)==20 and group.fixed_reference_pairs.nunique()==1
            denominator=int(group.fixed_reference_pairs.iloc[0]);assert denominator>0
            values=group.coassigned_pairs/denominator
            summaries.append(dict(anchor_left=left,anchor_right=right,within_reference_group=left==right,
                                  fixed_reference_pairs=denominator,mean_coassigned_fraction=float(values.mean()),
                                  min_coassigned_fraction=float(values.min()),max_coassigned_fraction=float(values.max()),
                                  mean_unresolved_fraction=float((group.unresolved_pairs/denominator).mean())))
        pd.DataFrame(summaries).to_csv(args.out/'contrast_summary.tsv',sep='\t',index=False)
        write_json(dict(original_seed_exact_labels=True,graph_unchanged=True,reference_pair_records_reconciled=4800,
                        half_stage_initial_equals_final=True,algorithmic_repeats=20,biological_specimens=1,
                        groups_min=min(x['supported_groups'] for x in fits),groups_max=max(x['supported_groups'] for x in fits),
                        no_new_CNA_truth=True,no_merged_or_RNA_fit=True,source=record,
                        estimand='Conditional coassignment of fixed reference cell pairs across Leiden seeds on one graph',
                        limitation='Algorithmic variation, not independent biological sampling; no support cutoff or confidence claim'),
                   args.out/'report.json')
        write_json({str(path.relative_to(args.out)):file_record(path) for path in sorted(args.out.rglob('*'))
                    if path.is_file() and path.name not in ('run_manifest.json','output_files.json')},args.out/'output_files.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
