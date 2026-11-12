"""Replay the exact executed DNA clustering block and record its merge decisions."""
import argparse
import ast
import copy
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import save_npz

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def assigned(node,name):
    return isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id==name for t in node.targets)


def replay_code(source):
    main=next(node for node in ast.parse(source).body if isinstance(node,ast.FunctionDef) and node.name=='main')
    loop=next(node for node in main.body if isinstance(node,ast.For) and isinstance(node.target,ast.Name) and node.target.id=='line')
    begin=next(i for i,node in enumerate(loop.body) if assigned(node,'a'))
    end=next(i for i,node in enumerate(loop.body) if assigned(node,'clone'))
    nodes=copy.deepcopy(loop.body[begin:end+1]);assert begin<end
    instrumented=[]
    for node in nodes:
        if isinstance(node,ast.While):
            body=[]
            for statement in node.body:
                if isinstance(statement,ast.Assign) and any(isinstance(t,ast.Subscript) and isinstance(t.value,ast.Name) and t.value.id=='lab' for t in statement.targets):
                    body.extend(ast.parse('_merge(dmin, i, j, lab)').body)
                body.append(statement)
                if assigned(statement,'pairs'):
                    body.extend(ast.parse('_pairs(pairs, lab)').body)
            node.body=body
        instrumented.append(node)
        if assigned(node,'lab'):
            instrumented.extend(ast.parse('_initial(lab)').body)
    return ast.fix_missing_locations(ast.Module(body=instrumented,type_ignores=[]))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    fit=args.task/'reconstruction_v1'
    manifest=json.loads((fit/'run_manifest.json').read_text());assert manifest['status']=='complete'
    record=manifest['input']['files']['DNA_source'];source=Path(record['path'])
    assert file_record(source)==record
    helper=Path(__file__).with_name('andor_dna_stability_audit.py')
    spec=importlib.util.spec_from_file_location('verified',helper)
    verified=importlib.util.module_from_spec(spec);spec.loader.exec_module(verified)
    verified.verified_run(fit)
    run=start_run(args.out,dict(script=Path(__file__),helper=helper,DNA_source=source,
                  fit=fit/'run_manifest.json',design=args.task/'CLUSTER_TRACE_DESIGN.md'),{},
                  dict(scope='Exact saved DNA clustering replay; no threshold selection or RNA inference'),mode='andor_dna_cluster_trace')
    try:
        spec=importlib.util.spec_from_file_location('executed_DNA',source)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        code=replay_code(source.read_text());(args.out/'executed_cluster_block.py').write_text(ast.unparse(code)+'\n')
        results={};stages={};graphs={};cell_orders={}
        for label in ('old','corrected'):
            directory=fit/label/'SNU-601';destination=args.out/label;destination.mkdir()
            with np.load(directory/'cn.npz',allow_pickle=True) as archive:
                z={key:archive[key] for key in ('cells','bins','cn','seg_ratio')}
            table=pd.read_csv(directory/'cells.tsv',sep='\t',dtype={'rg':str}).set_index('rg')
            assert len(set(z['cells']))==len(z['cells'])
            assert set(z['cells'])==set(table.index[table.pass_mapd])
            ann=pd.DataFrame(z['bins'],columns=['chromosome','start','end','arm'])
            ann[['start','end']]=ann[['start','end']].astype(np.int64)
            initial=[];rounds=[];merges=[]
            def capture_initial(lab):
                initial.append(lab.copy())
            def capture_pairs(pairs,lab):
                rounds.append(dict(round=len(rounds),group_sizes={str(k):int(np.sum(lab==k)) for k in np.unique(lab)},
                                   candidate_pairs=[dict(difference_Mb=float(d),left=int(i),right=int(j)) for d,i,j in pairs]))
            def capture_merge(d,i,j,lab):
                merges.append(dict(round=len(rounds)-1,difference_Mb=float(d),left=int(i),right=int(j),
                                   left_cells=int(np.sum(lab==i)),right_cells=int(np.sum(lab==j))))
            namespace=dict(module.__dict__,ad=ad,sc=sc,S=z['seg_ratio'],CN=z['cn'],ann=ann,
                           line='SNU-601',args=SimpleNamespace(min_clone=10,min_diff_mb=10.),
                           _initial=capture_initial,_pairs=capture_pairs,_merge=capture_merge)
            np.random.seed(module.SEED)
            exec(compile(code,str(source), 'exec'),namespace)
            assert len(initial)==1
            np.testing.assert_array_equal(namespace['clone'],table.loc[z['cells'],'clone'].to_numpy())
            frame=pd.DataFrame(dict(rg=z['cells'],initial_leiden=initial[0],final_group=namespace['clone']))
            frame.to_csv(destination/'stages.tsv',sep='\t',index=False)
            write_json(dict(rounds=rounds,merges=merges),destination/'merge_trace.json')
            a=namespace['a'];np.save(destination/'PCA.npy',a.obsm['X_pca'])
            save_npz(destination/'connectivities.npz',a.obsp['connectivities'])
            save_npz(destination/'distances.npz',a.obsp['distances'])
            graphs[label]=a.obsp['connectivities'].copy();cell_orders[label]=z['cells']
            results[label]=dict(exact_final_assignments_reproduced=True,cells=len(frame),
                                initial_groups=int(frame.initial_leiden.nunique()),
                                final_groups=int(frame.final_group.nunique()),merges=merges)
            stages[label]=frame.set_index('rg')
        assert set(stages['old'].index)==set(stages['corrected'].index)
        np.testing.assert_array_equal(cell_orders['old'],cell_orders['corrected'])
        difference=graphs['old']-graphs['corrected'];difference.eliminate_zeros()
        support_difference=(graphs['old']!=0)!=(graphs['corrected']!=0)
        joined=stages['old'].join(stages['corrected'],lsuffix='_old',rsuffix='_corrected')
        joined.to_csv(args.out/'joint_stages.tsv',sep='\t')
        pd.crosstab(joined.initial_leiden_old,joined.initial_leiden_corrected).to_csv(args.out/'initial_contingency.tsv',sep='\t')
        write_json(dict(replays=results,biological_specimens=1,
                        graph_comparison=dict(cells_aligned_exactly=True,
                            old_directed_nonzero_entries=int(graphs['old'].nnz),
                            corrected_directed_nonzero_entries=int(graphs['corrected'].nnz),
                            directed_entries_with_changed_support=int(support_difference.nnz),
                            directed_entries_with_changed_weight=int(difference.nnz),
                            maximum_absolute_weight_change=float(np.max(np.abs(difference.data))) if difference.nnz else 0.),
                        interpretation='Exact algorithm-stage attribution; neither grouping is certified biological truth'),args.out/'report.json')
        # Include nested stage arrays, sparse graphs and merge traces in the identity record.
        write_json({str(path.relative_to(args.out)):file_record(path) for path in sorted(args.out.rglob('*'))
                    if path.is_file() and path.name not in ('run_manifest.json','output_files.json')},args.out/'output_files.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
