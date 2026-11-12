"""Independent count/metric audit and descriptive collection of DNA perturbations."""
import argparse
from collections import Counter
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def verified_run(directory):
    manifest=json.loads((directory/'run_manifest.json').read_text())
    assert manifest['status']=='complete' and manifest['backend']=='cpu'
    for group in ('input','output'):
        for record in manifest[group]['files'].values():
            assert file_record(record['path'])['sha256']==record['sha256']
    inventory=directory/'output_files.json'
    if inventory.exists():
        for record in json.loads(inventory.read_text()).values():
            assert file_record(record['path'])['sha256']==record['sha256']
    return manifest


def independent_counts(path,cells,bins):
    """Chunked pandas/index assignment, independent of the streaming CSV reader."""
    ids=pd.Index(cells);keys=pd.MultiIndex.from_frame(bins[['chromosome','bin']])
    assert ids.is_unique and keys.is_unique
    matrix=np.zeros((len(cells),len(bins)),np.int64);seen=np.zeros_like(matrix,bool)
    totals=Counter();autosomes={f'chr{i}' for i in range(1,23)}
    for chunk in pd.read_csv(path,sep='\t',header=None,names=['rg','chromosome','bin','n'],
                             dtype={'rg':str,'chromosome':str,'bin':np.int64,'n':np.int64},chunksize=500000):
        rows=ids.get_indexer(chunk.rg);assert np.all(rows>=0)
        assert np.all(chunk.bin>=0) and np.all(chunk.n>0)
        for key,value in chunk.groupby(['rg','chromosome'],sort=False).n.sum().items():
            totals[key]+=int(value)
        keep=chunk.chromosome.isin(autosomes).to_numpy();selected=chunk.loc[keep]
        cols=keys.get_indexer(pd.MultiIndex.from_frame(selected[['chromosome','bin']]))
        rows=rows[keep];assert np.all(cols>=0)
        linear=rows*len(bins)+cols;assert len(np.unique(linear))==len(linear)
        assert not np.any(seen[rows,cols]);seen[rows,cols]=True
        matrix[rows,cols]=selected.n.to_numpy()
    return matrix,totals


def callable_flow(old,new,flow,bins,ann):
    keys={(row.chromosome,int(row.bin)):i for i,row in enumerate(bins.itertuples())}
    keep=np.zeros(len(bins),bool)
    for row in ann.itertuples():
        keep[keys[(row.chromosome,int(row.bin))]]=row.n_frac<=.1 and row.blacklist<=.1 and row.mappability>=.5
    entering=np.zeros(len(old),np.int64);leaving=entering.copy();retained=entering.copy();outside=entering.copy()
    for _,group in bins.groupby('chromosome',sort=False):
        columns=group.index.to_numpy();incoming=np.zeros(len(old),np.int64)
        # Reverse recurrence, independent of the producer's forward cumulative sum.
        for column in columns[::-1]:
            outgoing=old[:,column]+incoming-new[:,column]
            assert np.all(outgoing>=0) and np.all(outgoing<=old[:,column])
            np.testing.assert_array_equal(outgoing,flow[:,column]);incoming=outgoing
        assert np.all(incoming==0)
        source=keep[columns];target=np.r_[False,source[:-1]]
        for destination,mask in ((entering,~source&target),(leaving,source&~target),
                                 (retained,source&target),(outside,~source&~target)):
            destination+=flow[:,columns[mask]].sum(axis=1)
    np.testing.assert_array_equal(new[:,keep].sum(axis=1)-old[:,keep].sum(axis=1),entering-leaving)
    assert int(entering.sum()+leaving.sum()+retained.sum()+outside.sum())==int(flow.sum())
    return dict(entering_callable=entering,leaving_callable=leaving,
                staying_callable=retained,staying_outside_callable=outside)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    reference=Path(__file__).with_name('andor_dna_independent.py')
    spec=importlib.util.spec_from_file_location('independent_DNA',reference)
    metric=importlib.util.module_from_spec(spec);spec.loader.exec_module(metric)
    inputs=dict(script=Path(__file__),reference=reference,design=args.task/'DESIGN.md',
                audit_design=args.task/'AUDIT_DESIGN.md',preflight=args.task/'independent_preflight_v1/run_manifest.json')
    inputs['count_audit_preflight']=args.task/'audit_preflight_v1/run_manifest.json'
    for label in ('count_audit','reconstruction'):
        inputs[label]=args.task/(label+'_v1/run_manifest.json')
    for i in range(4):
        inputs[f'worker_{i}']=args.task/f'worker_{i}_v1/run_manifest.json'
    run=start_run(args.out,inputs,{},dict(scope='One exposed DNA specimen; independent conditional audit, not biological calibration'),
                  mode='andor_dna_stability_audit')
    try:
        verified_run(args.task/'independent_preflight_v1')
        verified_run(args.task/'audit_preflight_v1')
        count_dir=args.task/'count_audit_v1';counts_manifest=verified_run(count_dir)
        reconstruction=args.task/'reconstruction_v1';verified_run(reconstruction)
        with np.load(count_dir/'counts.npz') as archive:
            saved={key:archive[key] for key in archive.files}
        cells=saved['cells'];bins=pd.read_csv(count_dir/'bins.tsv',sep='\t')
        old,ot=independent_counts(Path(counts_manifest['input']['files']['old']['path']),cells,bins)
        new,nt=independent_counts(Path(counts_manifest['input']['files']['new']['path']),cells,bins)
        assert ot==nt
        np.testing.assert_array_equal(old,saved['old']);np.testing.assert_array_equal(new,saved['new'])
        original_annotation=Path(counts_manifest['input']['files']['annotation']['path'])
        assert (count_dir/'annotation.tsv').read_bytes()==original_annotation.read_bytes()
        ann=pd.read_csv(original_annotation,sep='\t')
        support=callable_flow(old,new,saved['flow'],bins,ann)
        moved=pd.DataFrame(dict(rg=cells,relocated_boundary_fragments=saved['flow'].sum(axis=1),**support))
        moved.to_csv(args.out/'coordinate_flow.tsv',sep='\t',index=False)
        comparison=json.loads((reconstruction/'comparison.json').read_text())
        checks=[metric.audit_comparison(reconstruction/'old/SNU-601',reconstruction/'corrected/SNU-601',comparison)]
        records=[];pairs=[];seen=set()
        for worker in range(4):
            directory=args.task/f'worker_{worker}_v1';manifest=verified_run(directory)
            report=json.loads((directory/'report.json').read_text())
            assert len(report['pairs'])==5
            for pair in report['pairs']:
                rep=pair['pair'];assert rep not in seen and 5*worker<=rep<5*(worker+1)
                seen.add(rep);assert pair['seed']==[20261108,1401,rep]
                folder=directory/f'pair_{rep:02d}'
                expected_A=np.random.default_rng(np.random.SeedSequence(pair['seed'])).binomial(new,.5)
                actual=[];fit_paths=[]
                for index,label in enumerate(('A','B')):
                    parent=folder/label;path=parent/'bins/SNU-601.tsv.gz'
                    assert file_record(path)==pair['half_counts'][index]
                    assert (parent/'bin_annotation_500kb.tsv').read_bytes()==original_annotation.read_bytes()
                    assert (parent/'bins/SNU-601.rg.txt').read_bytes()==(count_dir/'read_groups.txt').read_bytes()
                    matrix,_=independent_counts(path,cells,bins);actual.append(matrix)
                    np.testing.assert_array_equal(matrix,expected_A if label=='A' else new-expected_A)
                    fit_paths.append(parent/'DNA/SNU-601')
                np.testing.assert_array_equal(actual[0]+actual[1],new)
                assert pair['full_fragments']==int(new.sum())
                assert pair['A_fragments']==int(actual[0].sum()) and pair['B_fragments']==int(actual[1].sum())
                paths=[(reconstruction/'corrected/SNU-601',fit_paths[0],'A_to_full'),
                       (reconstruction/'corrected/SNU-601',fit_paths[1],'B_to_full'),
                       (fit_paths[0],fit_paths[1],'A_to_B')]
                full_denominator=pair['comparisons']['A_to_full']['anchor_cells']
                for left,right,name in paths:
                    result=pair['comparisons'][name]
                    checks.append(metric.audit_comparison(left,right,result))
                    records.append(dict(pair=rep,comparison=name,**result,
                        fixed_full_DNA_cells=full_denominator,
                        matched_assigned_fraction_of_fixed_full_DNA_cells=result['one_to_one_matched_assigned']/full_denominator))
                pairs.append(pair)
                print(f'Independent DNA pair audit {len(seen)}/20',flush=True)
        assert seen==set(range(20)) and len(checks)==61
        pd.DataFrame(records).to_csv(args.out/'comparison_metrics.tsv',sep='\t',index=False)
        reference_arms=sorted(set(ann.arm))
        support_rows=[]
        for arm in reference_arms:
            counts=Counter()
            for pair in pairs:
                flags=tuple(int(arm in pair['comparisons'][label+'_to_full']['candidate_inferred_difference_arms']) for label in ('A','B'))
                counts[flags]+=1
            support_rows.append(dict(arm=arm,pairs_00=counts[(0,0)],pairs_01=counts[(0,1)],
                                     pairs_10=counts[(1,0)],pairs_11=counts[(1,1)],
                                     supported_halves=counts[(0,1)]+counts[(1,0)]+2*counts[(1,1)],halves=40))
        pd.DataFrame(support_rows).to_csv(args.out/'arm_support.tsv',sep='\t',index=False)
        summaries=[]
        for name in ('A_to_full','B_to_full','A_to_B'):
            selected=[row for row in records if row['comparison']==name]
            values=[row['matched_assigned_fraction_of_fixed_full_DNA_cells'] for row in selected]
            summaries.append(dict(comparison=name,pairs=len(selected),
                                  mean_matched_assigned_fraction_of_fixed_full_DNA_cells=float(np.mean(values)),
                                  minimum_matched_fraction=min(values),maximum_matched_fraction=max(values),
                                  candidate_clone_counts=dict(Counter(row['candidate_clones'] for row in selected)),
                                  candidate_assigned_counts=[row['candidate_assigned'] for row in selected]))
        write_json(checks,args.out/'independent_checks.json')
        write_json(dict(comparisons_checked=len(checks),split_pairs=20,half_count_inputs_checked=40,
                        original_count_matrices_reconstructed=True,recursive_hashes_verified=True,
                        annotation_bytes_identical=True,autosomal_relocated_boundary_fragments=int(saved['flow'].sum()),
                        callable_flow_totals={key:int(value.sum()) for key,value in support.items()},
                        old_corrected_comparison=comparison,descriptive_summaries=summaries,
                        biological_specimens=1,independent_biological_validation=False,
                        interpretation='Conditional DNA measurement stability; no RNA membership or error-control claim'),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
