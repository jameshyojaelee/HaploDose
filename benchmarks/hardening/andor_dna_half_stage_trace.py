"""Attribute saved half-depth DNA coassignment to initial clustering or merging."""
import argparse
import ast
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def replay(module,code,directory,destination):
    """Saved S is exactly the float32 PCA input; saved CN is used without requantization."""
    with np.load(directory/'cn.npz',allow_pickle=True) as archive:
        z={name:archive[name] for name in ('cells','bins','cn','seg_ratio')}
    table=pd.read_csv(directory/'cells.tsv',sep='\t',dtype={'rg':str,'clone':str},keep_default_na=False).set_index('rg')
    assert table.index.is_unique and pd.Index(z['cells']).is_unique
    assert set(z['cells'])==set(table.index[table.pass_mapd])
    ann=pd.DataFrame(z['bins'],columns=['chromosome','start','end','arm'])
    ann[['start','end']]=ann[['start','end']].astype(np.int64)
    initial=[];rounds=[];merges=[]
    def capture_pairs(pairs,lab):
        rounds.append(dict(round=len(rounds),candidate_pairs=[dict(difference_Mb=float(d),left=int(i),right=int(j)) for d,i,j in pairs]))
    def capture_merge(d,i,j,lab):
        merges.append(dict(round=len(rounds)-1,difference_Mb=float(d),left=int(i),right=int(j),
                           left_cells=int(np.sum(lab==i)),right_cells=int(np.sum(lab==j))))
    namespace=dict(module.__dict__,ad=ad,sc=sc,S=z['seg_ratio'],CN=z['cn'],ann=ann,
                   line='SNU-601',args=SimpleNamespace(min_clone=10,min_diff_mb=10.),
                   _initial=lambda lab:initial.append(lab.copy()),_pairs=capture_pairs,_merge=capture_merge)
    np.random.seed(module.SEED)
    exec(compile(code,module.__file__,'exec'),namespace)
    assert len(initial)==1
    np.testing.assert_array_equal(namespace['clone'],table.loc[z['cells'],'clone'].to_numpy())
    # Apply the executed minimum-group-size rule to both stages. Small initial
    # groups remain unresolved; they must never count as an assigned clone.
    size=pd.Series(initial[0]).value_counts()
    before=np.array([f'initial_{k}' if size[k]>=10 else 'unassigned' for k in initial[0]],object)
    frame=pd.DataFrame(dict(rg=z['cells'],initial_leiden=initial[0],initial_group=before,
                            final_group=namespace['clone'])).set_index('rg')
    result=pd.DataFrame(index=table.index)
    for column in ('initial_group','final_group'):
        result[column]=frame[column].reindex(table.index,fill_value='excluded')
    np.testing.assert_array_equal(result.final_group,table.clone)
    destination.mkdir()
    frame.to_csv(destination/'stages.tsv',sep='\t')
    write_json(dict(rounds=rounds,merges=merges),destination/'merge_trace.json')
    return result,dict(exact_final_assignments_reproduced=True,profile_cells=len(frame),
                       initial_leiden_groups=int(len(size)),initial_supported_groups=int((size>=10).sum()),
                       final_groups=int(table.loc[~table.clone.isin(['excluded','unassigned']),'clone'].nunique()),
                       accepted_merges=len(merges))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    root=Path(__file__).parent
    trace=load('exact_stage_trace',root/'andor_dna_cluster_trace.py')
    contrast=load('fixed_pair_counts',root/'andor_dna_contrast_stability.py')
    verified=load('verified_stage_inputs',root/'andor_dna_stability_audit.py')
    inputs=dict(script=Path(__file__),trace_helper=root/'andor_dna_cluster_trace.py',
                contrast_helper=root/'andor_dna_contrast_stability.py',
                verification_helper=root/'andor_dna_stability_audit.py',
                design=args.task/'HALF_STAGE_DESIGN.md')
    for name in ('reconstruction_v1','audit_v1','contrast_preflight_v1','cluster_trace_v2'):
        inputs[name]=args.task/name/'run_manifest.json'
    for worker in range(4):
        inputs[f'worker_{worker}']=args.task/f'worker_{worker}_v1/run_manifest.json'
    run=start_run(args.out,inputs,{},dict(scope='Conditional stage attribution in one exposed DNA specimen; no RNA inference'),
                  mode='andor_dna_half_stage_trace')
    try:
        for name in ('audit_v1','contrast_preflight_v1','cluster_trace_v2'):
            verified.verified_run(args.task/name)
        full=verified.verified_run(args.task/'reconstruction_v1')
        # First replay the two already traced fits as a gate on this new wrapper.
        record=full['input']['files']['DNA_source'];source=Path(record['path'])
        assert file_record(source)==record
        module=load('executed_full_DNA',source);code=trace.replay_code(source.read_text())
        gates={}
        for label in ('old','corrected'):
            _,gates[label]=replay(module,code,args.task/f'reconstruction_v1/{label}/SNU-601',args.out/f'gate_{label}')
        anchor=contrast.read_groups(args.task/'reconstruction_v1/corrected/SNU-601/cells.tsv')
        fits=[];pairs=[];sources={}
        for worker in range(4):
            directory=args.task/f'worker_{worker}_v1';manifest=verified.verified_run(directory)
            record=manifest['input']['files']['DNA_source'];source=Path(record['path'])
            assert file_record(source)==record
            sources[str(worker)]=record
            module=load(f'executed_half_DNA_{worker}',source);code=trace.replay_code(source.read_text())
            (args.out/f'executed_cluster_block_{worker}.py').write_text(ast.unparse(code)+'\n')
            for rep in range(worker*5,(worker+1)*5):
                for half in ('A','B'):
                    frame,metrics=replay(module,code,directory/f'pair_{rep:02d}/{half}/DNA/SNU-601',
                                         args.out/f'pair_{rep:02d}_{half}')
                    assert set(frame.index)==set(anchor.index)
                    frame=frame.loc[anchor.index]
                    metrics.update(pair=rep,half=half);fits.append(metrics)
                    for stage,column in (('initial','initial_group'),('final','final_group')):
                        records=contrast.pair_counts(anchor.tolist(),frame[column].tolist())
                        for item in records:
                            item.update(pair=rep,half=half,stage=stage)
                        pairs.extend(records)
                    print(json.dumps(metrics),flush=True)
        table=pd.DataFrame(pairs)
        table.to_csv(args.out/'stage_pair_counts.tsv',sep='\t',index=False)
        pd.DataFrame(fits).to_csv(args.out/'fit_stages.tsv',sep='\t',index=False)
        index=['pair','half','anchor_left','anchor_right']
        before=table.loc[table.stage=='initial'].set_index(index)
        after=table.loc[table.stage=='final'].set_index(index).loc[before.index]
        for column in ('fixed_reference_pairs','unresolved_pairs','both_assigned_pairs'):
            np.testing.assert_array_equal(before[column],after[column])
        assert (after.coassigned_pairs>=before.coassigned_pairs).all()
        delta=after[['fixed_reference_pairs','within_reference_group']].copy()
        delta['merge_added_coassigned_pairs']=after.coassigned_pairs-before.coassigned_pairs
        delta.to_csv(args.out/'merge_pair_changes.tsv',sep='\t')
        write_json(dict(gates=gates,half_fits=len(fits),all_final_assignments_exact=True,
                        fits_with_accepted_merges=sum(fit['accepted_merges']>0 for fit in fits),
                        total_accepted_merges=sum(fit['accepted_merges'] for fit in fits),
                        biological_specimens=1,independent_pair_seeds=20,complementary_halves_dependent=True,
                        source_records=sources,reference_groups='Uncertain full-depth DNA reconstruction, same molecules',
                        estimand='Stage-specific coassignment of fixed full-depth reference cell pairs',
                        limitation='No certified clones, RNA accuracy, biological error calibration or selected support cutoff'),
                   args.out/'report.json')
        write_json({str(path.relative_to(args.out)):file_record(path) for path in sorted(args.out.rglob('*'))
                    if path.is_file() and path.name not in ('run_manifest.json','output_files.json')},args.out/'output_files.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
