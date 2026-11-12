"""Independent recovery recount, immutable outputs, and membership-stage diagnosis."""
import argparse
import itertools
import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record


def verify(record):
    if file_record(record['path'])['sha256']!=record['sha256']:raise AssertionError('Recorded artifact changed')


def recovery(cells,diff):
    s1=cells[cells.specimen=='S1'];normal=s1.truth_clone=='normal'
    normal_fpr=float(s1.loc[normal,'malignant_call'].sum()/normal.sum())
    tumour_tpr=float(s1.loc[~normal,'malignant_call'].sum()/(~normal).sum())
    sub=s1[s1.truth_clone.isin(['clone1','clone2']) & s1.malignant_call]
    # Reproduce the original unusual two-row diagonal estimand explicitly; its
    # separate exactly-two-groups check must pass too. No column omission.
    rows=sorted(sub.truth_clone.unique());cols=sorted(sub.clone.unique())
    table=np.array([[int(((sub.truth_clone==r)&(sub.clone==c)).sum()) for c in cols] for r in rows])
    diagonal=max(sum(table[i,i] for i in range(min(table.shape))),
                 sum(table[len(rows)-1-i,i] for i in range(min(table.shape)))) if table.size else 0
    accuracy=diagonal/len(sub) if len(sub) else None
    use=(diff.specimen=='S1') & ~diff.clone_a.str.endswith('_c0') & ~diff.clone_b.str.endswith('_c0')
    values=dict(normal_fpr=normal_fpr,tumour_tpr=tumour_tpr,
        s2_one_group=len(set(cells.loc[(cells.specimen=='S2')&cells.malignant_call,'clone'])),
        s1_two_groups=len(cols),assignment_accuracy=accuracy,difference_chromosomes=sorted(set(diff.loc[use,'chromosome'])))
    passes=dict(normal_fpr=normal_fpr<.1,tumour_tpr=tumour_tpr>.9,s2_one_group=values['s2_one_group']==1,
        s1_two_groups=len(cols)==2,assignment_accuracy=accuracy is not None and accuracy>.65,
        difference_chromosomes=values['difference_chromosomes']==['chr12'])
    return values,{key:bool(value) for key,value in passes.items()}


def membership(truth,labels):
    truth=np.asarray(truth);labels=np.asarray(labels);use=np.isin(truth,['clone1','clone2'])
    truth,labels=truth[use],labels[use];observed=[v for v in np.unique(labels) if v>=0]
    # Dummy labels permit one group; unmatched/extra groups and abstention count
    # as incorrect. Labels are aligned for scoring only, never sent to inference.
    choices=observed+[-2,-3];best=0
    for a,b in itertools.permutations(choices,2):
        best=max(best,int(((truth=='clone1')&(labels==a)).sum()+((truth=='clone2')&(labels==b)).sum()))
    return dict(n=len(truth),one_to_one_accuracy=best/len(truth) if len(truth) else None,
        groups=len(observed),assigned_fraction=float((labels>=0).mean()) if len(truth) else None)


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--attempt',default='v1',choices=['v1','v2'])
    a=p.parse_args();a.out.mkdir(exist_ok=False);inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md')
    for arm in ('hard','margin'):
        inputs[arm]=a.task/f'{arm}_{a.attempt}/run_manifest.json'
        for seed in (20261031,20261107,20261108,20261109,20261110):
            case=a.task/f'{arm}_{a.attempt}/seed_{seed}'
            for name in ('trace.json','summary.json','scored_cells.tsv','fit/run_manifest.json'):
                inputs[f'{arm}:{seed}:{name}']=case/name
    rec=start_run(a.out,inputs,dict(arms=['hard','margin']),dict(scope='Independent recount; no biological claim'),mode='joint_recovery_audit')
    try:
        rows=[];stages=[]
        for arm in ('hard','margin'):
            directory=a.task/f'{arm}_{a.attempt}';manifest=json.loads(inputs[arm].read_text())
            if manifest['status']!='complete':raise AssertionError('Incomplete arm')
            for r in manifest['input']['files'].values():verify(r)
            for r in manifest['output']['files'].values():verify(r)
            summaries=json.loads((directory/'summary.json').read_text())['cases']
            if [r['seed'] for r in summaries]!=[20261031,20261107,20261108,20261109,20261110]:raise AssertionError('Cases changed')
            for summary in summaries:
                seed=summary['seed'];case=directory/f'seed_{seed}';fit=case/'fit'
                fm=json.loads((fit/'run_manifest.json').read_text())
                if fm['status']!='complete':raise AssertionError('Incomplete fit')
                for r in fm['input']['files'].values():verify(r)
                for r in fm['output']['files'].values():verify(r)
                qc=json.loads((fit/'qc.json').read_text())
                if qc['confirmation_status']!='joint_not_independent':raise AssertionError('False confirmation claim')
                cells=pd.read_csv(case/'scored_cells.tsv',sep='\t');diff=pd.read_csv(fit/'differences.tsv',sep='\t')
                values,passes=recovery(cells,diff)
                for key,value in values.items():
                    previous=summary['criteria'][key]
                    if previous['value']!=value or previous['passed']!=passes[key]:raise AssertionError('Criterion recount failed')
                if summary['all_numeric_criteria_pass']!=all(passes.values()):raise AssertionError('Wrong success status')
                rows.append(dict(arm=arm,seed=seed,values=values,passes=passes))
                source=Path(manifest['input']['files'][f'{seed}:counts.h5']['path'])
                with h5py.File(source,'r') as h:ids=h['cell_ids'].asstr()[:]
                truth_path=Path(manifest['input']['files'][f'{seed}:cells.tsv']['path'])
                truth=pd.read_csv(truth_path,sep='\t').set_index('cell_id').loc[ids]
                exported=pd.read_csv(fit/'cells.tsv.gz',sep='\t')
                for column in ('cell_id','specimen','clone','malignant_call'):
                    missing=cells[column].isna()
                    np.testing.assert_array_equal(missing,exported[column].isna())
                    np.testing.assert_array_equal(cells.loc[~missing,column],exported.loc[~missing,column])
                np.testing.assert_array_equal(cells.truth_clone,truth.loc[cells.cell_id,'truth_clone'])
                for index,event in enumerate(json.loads((case/'trace.json').read_text())):
                    if event['stage'] not in {'joint_proposal','merge_membership'}:continue
                    specimen=truth.iloc[event['specimen_rows']]
                    if set(specimen.specimen)!= {'S1'}:continue
                    if event['stage']=='joint_proposal':
                        if 'labels' not in event:continue
                        selected=specimen.iloc[event['proposal_cells']].truth_clone.to_numpy()
                        stages.append(dict(arm=arm,seed=seed,index=index,stage='joint_proposal',source=event['source'],
                            original=membership(selected,event['original_labels']),joint=membership(selected,event['labels']),
                            mean_joint_margin=float(np.abs(2*np.asarray(event['membership'])-1).mean())))
                    else:
                        result=dict(arm=arm,seed=seed,index=index,stage='merge_membership',pairs=event['pairs'])
                        for field in ('full_responsibility','proposal_responsibility'):
                            R=np.asarray(event[field]);labels=R.argmax(axis=1);labels[labels==0]=-1
                            result[field]=membership(specimen.truth_clone.to_numpy(),labels)
                            result[field]['mean_max_probability']=float(R.max(axis=1).mean())
                        stages.append(result)
        baseline=json.loads((a.task/'baseline/source_manifest.json').read_text());core=0
        for name,digest in baseline.items():
            if file_record(a.task/'baseline'/name)['sha256']!=digest:raise AssertionError('Baseline changed')
            if name.startswith(('python/','tests/','crates/')):
                if file_record(a.task.parents[1]/name)['sha256']!=digest:raise AssertionError('Core changed')
                core+=1
        report=dict(cases=rows,membership_stages=stages,baseline_files=len(baseline),core_files=core,
            decision='No calibration or original-thinning success claim; retain all failures and assess complete screen before further sweep')
        (a.out/'report.json').write_text(json.dumps(report,indent=2));finish_run(a.out,rec,backend='cpu')
        print(json.dumps(dict(cases=len(rows),membership_records=len(stages),baseline_files=len(baseline),core_files=core)),flush=True)
    except Exception as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
