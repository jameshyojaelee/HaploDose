"""Check whether label-corrected numerical changes reach accepted RNA groups."""
import argparse
import json
from pathlib import Path

import pandas as pd

from haplodose.provenance import start_run, finish_run, file_record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    args=p.parse_args();args.out.mkdir(exist_ok=False)
    prior=args.task.parents[1]/'dev/rna_population_baseline_20261109'
    inputs=dict(script=Path(__file__),design=args.task/'LABEL_DESIGN.md',
        labels=args.task/'label_v1/all_cases.tsv',label_report=args.task/'label_v1/report.json',
        solvers=args.task/'probe_v3/solver_runs.tsv',probe=args.task/'probe_v3/report.json',
        decisions=prior/'summary_v2/decisions.tsv',original_report=prior/'summary_v2/report.json')
    for name,directory in (('labels',args.task/'label_v1'),('solvers',args.task/'probe_v3'),('original',prior/'summary_v2')):
        path=directory/'run_manifest.json';manifest=json.loads(path.read_text());assert manifest['status']=='complete'
        for record in manifest['output']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        inputs[name+':manifest']=path
    run=start_run(args.out,inputs,dict(new_fits=0,new_tests=0),
        dict(scope='Consequence audit of unchanged decisions; no new calibration or biological claim'),mode='rna_population_optimizer_summary')
    try:
        labels=pd.read_csv(inputs['labels'],sep='\t');decisions=pd.read_csv(inputs['decisions'],sep='\t')
        assert len(labels)==998 and len(decisions)==600
        assert not labels[['role','replicate']].duplicated().any()
        assert not decisions[['role','replicate']].duplicated().any()
        joint=decisions.merge(labels,on=['role','replicate'],validate='one_to_one',how='left',indicator=True)
        assert joint['_merge'].eq('both').all()
        assert set(joint.role)=={'eval_nb','eval_state','eval_two'}
        original=json.loads(inputs['original_report'].read_text())
        rows=[]
        for role,frame in joint.groupby('role'):
            assert len(frame)==200
            reported=next(r for r in original['endpoints'] if r['role']==role)
            accepted=frame.accepted
            assert int(accepted.sum())==reported['events']
            changed=frame.label_invariant_membership_changes.gt(0)
            rows.append(dict(role=role,specimens=200,accepted=int(accepted.sum()),
                changed_partitions=int(changed.sum()),accepted_with_changed_partition=int((accepted&changed).sum()),
                accepted_with_literal_label_change=int((accepted&frame.literal_label_changes.gt(0)).sum()),
                accepted_with_pure_label_reversal=int((accepted&frame.pure_label_reversal).sum())))
        solvers=pd.read_csv(inputs['solvers'],sep='\t');assert len(solvers)==32
        assert solvers[['role','replicate','start','solver']].drop_duplicates().shape[0]==32
        report=dict(evaluation=rows,scientific_event_counts_exact=True,new_fits=0,new_tests=0,
            corrected_probe_cases=len(solvers[['role','replicate']].drop_duplicates()),
            corrected_probe_converged_runs=int(solvers.success.sum()),
            corrected_probe_material_supported_gains=int(solvers.material_supported_gain_over_saved_fine.sum()),
            historical_literal_change_cases=67,
            label_invariant_change_cases=int(labels.label_invariant_membership_changes.gt(0).sum()),
            pure_reversal_cases=int(labels.pure_label_reversal.sum()),
            numerical_global_optimality_certified=False,scientific_criteria_changed=False,production_changed=False,
            biological_validation=False,method_promoted=False)
        pd.DataFrame(rows).to_csv(args.out/'decision_consequences.tsv',sep='\t',index=False)
        (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        for record in run['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(args.out,run,backend='cpu')
    except Exception as error:finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
