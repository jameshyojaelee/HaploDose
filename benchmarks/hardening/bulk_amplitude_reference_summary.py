"""Audit-gated, failure-inclusive reference-sensitivity summary; no promotion."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md')
    for i in range(9):
        inputs[f'panel_{i}']=a.task/f'refit_{i:03d}_v1/report.json'
        inputs[f'audit_{i}']=a.task/f'audit_{i:03d}_v1/report.json'
        inputs[f'effects_{i}']=a.task/f'audit_{i:03d}_v1/paired_effects.json'
        inputs[f'scores_{i}']=a.task/f'refit_{i:03d}_v1/patient_scores.tsv'
        for role in ('panel','audit'):
            directory=a.task/(f'refit_{i:03d}_v1' if role=='panel' else f'audit_{i:03d}_v1')
            manifest=json.loads((directory/'run_manifest.json').read_text());assert manifest['status']=='complete'
            for record in manifest['output']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
    rec=start_run(a.out,inputs,{},dict(scope='Descriptive exposed-panel sensitivity,not independent validation or calibrated uncertainty'),mode='bulk_amplitude_reference_summary')
    try:
        rows=[];marginal=[]
        for i in range(9):
            report=json.loads(inputs[f'panel_{i}'].read_text());audit=json.loads(inputs[f'audit_{i}'].read_text())
            assert report['reference_refit']==i and audit['all_checks_passed']
            assert audit['normal_source_correspondence_checks']==20
            effects=json.loads(inputs[f'effects_{i}'].read_text());assert len(effects)==80
            scores=pd.read_csv(inputs[f'scores_{i}'],sep='\t',float_precision='round_trip')
            for row in effects:
                assert row['total_comparison_patients']==24
                assert len(row['patients'])==row['paired_patients']
                # Independently recount the paired mean from the patient table.
                one=scores[(scores.method==row['method'])&(scores.split=='evaluation')].set_index('sample')
                two=scores[(scores.method==row['control'])&(scores.split=='evaluation')].set_index('sample')
                endpoint=row['endpoint'];ids=one.index[one[endpoint].notna()&two[endpoint].notna()]
                assert list(ids)==row['patients']
                value=sum(float(one.at[s,endpoint]-two.at[s,endpoint]) for s in ids)/len(ids) if len(ids) else np.nan
                recorded=row['mean_paired_effect']
                guard=len(ids)*np.finfo(float).eps*(1+abs(value))
                assert (not len(ids) and recorded is None) or abs(value-recorded)<=guard
                rows.append(dict(reference_refit=i,**row))
            for method,group in scores[scores.split=='evaluation'].groupby('method'):
                for endpoint in ('pearson','recall_loss','fpr_loss','precision_loss'):
                    marginal.append(dict(reference_refit=i,method=method,endpoint=endpoint,
                        endpoint_patients=int(group[endpoint].notna().sum()),total_patients=len(group),mean=group[endpoint].mean()))
        frame=pd.DataFrame(rows);frame.drop(columns=['patients']).to_csv(a.out/'paired_effects_by_reference.tsv',sep='\t',index=False)
        pd.DataFrame(marginal).to_csv(a.out/'marginal_means_by_reference.tsv',sep='\t',index=False)
        stability=[]
        for (method,control,endpoint),group in frame.groupby(['method','control','endpoint']):
            assert len(group)==9
            boot=group[group.reference_refit>0];direction=-1 if endpoint=='fpr_loss' else 1
            values=direction*boot.mean_paired_effect.to_numpy(float)
            stability.append(dict(method=method,control=control,endpoint=endpoint,
                original_effect=float(group.loc[group.reference_refit==0,'mean_paired_effect'].iloc[0]),
                bootstrap_reference_draws=8,finite_draws=int(np.isfinite(values).sum()),
                favourable_draws=int((values>0).sum()),unfavourable_draws=int((values<0).sum()),
                zero_draws=int((values==0).sum()),minimum_effect=float(boot.mean_paired_effect.min()),
                maximum_effect=float(boot.mean_paired_effect.max()),minimum_paired_patients=int(boot.paired_patients.min()),
                maximum_paired_patients=int(boot.paired_patients.max())))
        pd.DataFrame(stability).to_csv(a.out/'reference_stability.tsv',sep='\t',index=False)
        (a.out/'report.json').write_text(json.dumps(dict(reference_fits=9,bootstrap_reference_draws=8,
            independent_reference_cohorts=1,comparison_patients=24,calibration_patients=24,
            fixed_candidates=2,fixed_comparators=10,endpoints=4,paired_contrasts_per_reference=80,
            new_hypothesis_tests=0,confidence_intervals='not estimated at reference-sensitivity stage',
            precision_estimand='Conditional on both methods making classified calls; paired populations can change with reference',
            all_nine_independent_audits_passed=True,method_promoted=False),indent=2)+'\n')
        finish_run(a.out,rec,backend='CPU scalar paired summaries')
    except BaseException as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
