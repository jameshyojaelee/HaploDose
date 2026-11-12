"""Audit the fixed simple RNA baseline, rank decisions and joint error/power."""
import argparse
from fractions import Fraction
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import beta

from haplodose.clones import CloneParams
from haplodose.provenance import start_run, finish_run, file_record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    args=p.parse_args();args.out.mkdir(exist_ok=False)
    source=Path(__file__).parent
    spec=importlib.util.spec_from_file_location('independent_exact_population_bounds',source/'joint_fresh_summary.py')
    numerical=importlib.util.module_from_spec(spec);spec.loader.exec_module(numerical);numerical.FAMILY=156
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',exact_bounds=source/'joint_fresh_summary.py',
        preflight=args.task/'preflight_v2/checks.json')
    history=args.task.parents[1]/'dev/joint_fresh_null_20261108/summary_v1'
    history_manifest=json.loads((history/'run_manifest.json').read_text())
    assert history_manifest['status']=='complete'
    history_record=history_manifest['output']['files']['summary.json']
    assert file_record(history_record['path'])['sha256']==history_record['sha256']
    inputs['historical_summary']=Path(history_record['path'])
    inputs['historical_manifest']=history/'run_manifest.json'
    rows={}
    for part in ('calibration','evaluation'):
        directory=args.task/(part+'_v2');m=json.loads((directory/'run_manifest.json').read_text())
        assert m['status']=='complete'
        for record in m['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        for record in m['output']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        inputs[part+':manifest']=directory/'run_manifest.json'
        inputs[part+':records']=directory/'records.jsonl'
        for line in (directory/'records.jsonl').read_text().splitlines():
            row=json.loads(line);rows.setdefault(row['role'],[]).append(row)
    run=start_run(args.out,inputs,dict(confidence_family=156,alpha=.05,threshold_selection=False),
        dict(scope='Conditional fixed-bank simulated error and power; no biological calibration'),mode='rna_population_summary')
    try:
        assert all(json.loads(inputs['preflight'].read_text()).values())
        assert CloneParams().min_cells==30
        assert {k:len(v) for k,v in rows.items()}==dict(cal_nb=199,cal_state=199,eval_nb=200,eval_state=200,eval_two=200)
        for role,values in rows.items():
            assert [v['replicate'] for v in values]==list(range(len(values)))
            for row in values:
                assert file_record(row['arrays']['path'])['sha256']==row['arrays']['sha256']
                with np.load(row['arrays']['path']) as a:
                    probability=a['probability'];truth=a['simulated_CNA_labels'];labels=probability.argmax(axis=1)
                    assert probability.shape==(600,2) and truth.shape==(600,)
                    assert np.isfinite(probability).all() and (probability>=0).all()
                    # Softmax unit-sum floating error derives from binary64 operations.
                    assert np.max(np.abs(probability.sum(axis=1)-1)) <= 32*np.finfo(float).eps
                    assert np.bincount(labels,minlength=2).min()==row['minimum_MAP_group']
                    if role=='eval_two':
                        table=[[int(np.sum((truth==t)&(labels==k))) for k in (0,1)] for t in (0,1)]
                        matching=max(sum(table[t][assignment[t]] for t in (0,1)) for assignment in itertools.permutations((0,1)))/600
                        assert matching==row['all_cell_two_group_correctness']
                for fit in row['fits']:assert fit['independent_loglik_difference']<=fit['roundoff_budget']
        def statistic(row,key):return row[key] if row['minimum_MAP_group']>=30 else 0.
        def decisions(observed_key,calibration_key):
            banks={role:[statistic(r,calibration_key) for r in rows[role]] for role in ('cal_nb','cal_state')}
            result={}
            for role in ('eval_nb','eval_state','eval_two'):
                result[role]=[]
                for row in rows[role]:
                    score=statistic(row,observed_key)
                    fractions=[Fraction(1+sum(x>=score for x in bank),200) for bank in banks.values()]
                    rank=max(fractions);accepted=rank<=Fraction(1,20) and row['minimum_MAP_group']>=30
                    # Independent strict order-statistic formulation, including tied scores.
                    expected=score>max(sorted(bank,reverse=True)[9] for bank in banks.values()) and row['minimum_MAP_group']>=30
                    assert accepted==expected
                    result[role].append(dict(role=role,replicate=row['replicate'],score=score,
                        accepted=accepted,p_numerator=rank.numerator,p_denominator=rank.denominator))
            return result,banks
        primary,banks=decisions('score','score')
        refinements=[decisions(o,c)[0] for o,c in itertools.product(('score','coarse_score'),repeat=2)]
        differences={role:sum(any(other[role][i]['accepted']!=row['accepted'] for other in refinements)
            for i,row in enumerate(primary[role])) for role in primary}
        endpoints=[]
        for role in ('eval_nb','eval_state','eval_two'):
            k=sum(row['accepted'] for row in primary[role]);power=role=='eval_two'
            independent=numerical.bound(k,200,power)
            scipy=(0. if k==0 else float(beta.ppf(.05/156,k,201-k))) if power else (1. if k==200 else float(beta.isf(.05/156,k+1,200-k)))
            # Compare by forward tail residual and report discrepancies rather than
            # inventing a universal numerical epsilon for a scientific boundary.
            endpoints.append(dict(role=role,events=k,n=200,rate=k/200,
                simultaneous_lower=independent if power else None,
                simultaneous_upper=None if power else independent,
                scipy_reference=scipy,independent_reference_discrepancy=abs(independent-scipy),
                null_error_requirement_met=None if power else bool(independent<=.05)))
        historical=json.loads(inputs['historical_summary'].read_text())
        old_rows=historical['summary']+historical['retained_historical']
        assert historical['family']==153 and len(old_rows)==153
        retained=[];cache={}
        for row in old_rows:
            key=(row['events'],row['n'],row['direction']=='lower')
            if key not in cache:cache[key]=numerical.bound(*key)
            retained.append(dict(**{k:v for k,v in row.items() if k not in
                ('family','simultaneous_bound','previous_bound','error_criterion_met')},
                family=156,simultaneous_bound=cache[key],
                previous_bound=dict(family=153,value=row['simultaneous_bound']),
                error_criterion_met=None if key[2] else cache[key]<=.05))
        assert len(retained)+len(endpoints)==156
        (args.out/'retained_historical_bounds.json').write_text(json.dumps(retained,indent=2)+'\n')
        alternative=rows['eval_two'];accepted=np.array([r['accepted'] for r in primary['eval_two']])
        accuracy=np.array([r['all_cell_two_group_correctness'] for r in alternative])
        report=dict(endpoints=endpoints,confidence_family=156,confidence='Simultaneous one-sided95% exact; fixed realized banks',
            all_cell_uncalibrated_two_group_correctness=float(accuracy.mean()),
            accepted_two_group_cell_coverage=float(accepted.mean()),
            correct_assignments_over_all_cells_retaining_abstention=float(np.mean(accuracy*accepted)),
            correctness_conditional_on_acceptance=float(accuracy[accepted].mean()) if accepted.any() else None,
            rank_decision_disagreements_under_refinement=differences,
            maximum_score_refinement_change=max(r['refinement_score_change'] for v in rows.values() for r in v),
            cases_with_MAP_refinement_change=sum(r['refinement_MAP_changes']>0 for v in rows.values() for r in v),
            minimum_MAP_group_support=30,calibration_bank_size_each=199,
            statistic='Complete PCA plus fitted tied-variance mixture likelihood gain; support veto applied in both observed and null scores',
            oracle_inputs_used=False,all_cells_hard_assigned_before_specimen_gate=True,
            no_reference_like_or_normal_cell_detection=True,biological_specimens=0,
            unsupported_genomic_clone_claim_withheld=True,method_promoted=False,
            historical_153_family_and_original_recovery_criteria_unchanged=True,
            historical_endpoints_recomputed_in_new_156_family=len(retained))
        pd.DataFrame([r for v in primary.values() for r in v]).to_csv(args.out/'decisions.tsv',sep='\t',index=False)
        pd.DataFrame(endpoints).to_csv(args.out/'endpoints.tsv',sep='\t',index=False)
        (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        for record in run['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
