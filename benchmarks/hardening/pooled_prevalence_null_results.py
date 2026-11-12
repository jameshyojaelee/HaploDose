"""Read the prespecified n=10 null look after both complete null panels finish."""
import argparse
from fractions import Fraction
import importlib.util
import json
from pathlib import Path

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--pilot',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    helper=Path(__file__).with_name('pooled_prevalence_pilot_summary.py')
    spec=importlib.util.spec_from_file_location('frozen_null_summary',helper)
    exact=importlib.util.module_from_spec(spec);spec.loader.exec_module(exact)
    assert exact.ALPHA==Fraction(1,80) and exact.TARGET==Fraction(1,20)
    inputs=dict(script=Path(__file__),helper=helper,design=args.task/'DESIGN.md',
        original_design=args.pilot/'PILOT_DESIGN.md',
        confidence_numerics=args.pilot/'numerics_v1/report.json')
    for condition in('matched','state'):
        for start in(0,5):
            folder=args.pilot/f'pilot_{condition}_{start:02d}_v1'
            inputs[f'{condition}{start}_results']=folder/'results.json'
            inputs[f'{condition}{start}_manifest']=folder/'run_manifest.json'
    rec=start_run(args.out,inputs,dict(phase=8,n=10,primary_nulls=2,
        total_planned_looks=2,one_sided_alpha=str(exact.ALPHA),
        target=str(exact.TARGET),planned_final_n=100),
        dict(scope='Same planned n10 null look; power pending, biological n=0'),
        mode='pooled_prevalence_prespecified_null_look')
    try:
        assert json.loads(inputs['confidence_numerics'].read_text())['exact_rational_confidence_inversion']
        checks=exact.checks()
        results=[]
        for condition in('matched','state'):
            for start in(0,5):
                assert json.loads(inputs[f'{condition}{start}_manifest'].read_text())['status']=='complete'
                batch=json.loads(inputs[f'{condition}{start}_results'].read_text())
                assert len(batch)==10
                assert all(r['condition']==condition and r['phase']==8 and
                    r['assignment_n']==600 and isinstance(r['null_error'],bool)for r in batch)
                for arm in('production','pooled_prevalence'):
                    assert sorted(r['replicate']for r in batch if r['arm']==arm)==list(range(start,start+5))
                results.extend(batch)
        assert len(results)==40
        summaries=[]
        for condition in('matched','state'):
            arms={}
            for arm in('production','pooled_prevalence'):
                rows=sorted((r for r in results if r['condition']==condition and r['arm']==arm),
                    key=lambda r:r['replicate'])
                assert [r['replicate']for r in rows]==list(range(10))
                assert len({r['seed']for r in rows})==10
                arms[arm]=rows
                errors=sum(r['null_error']for r in rows)
                bound=exact.upper(errors,10,exact.ALPHA)
                excess=1-exact.cdf(errors-1,10,exact.TARGET)if errors else Fraction(1)
                summaries.append(dict(condition=condition,arm=arm,workflows=10,
                    false_splits=errors,simultaneous_false_split_upper=bound['value'],
                    simultaneous_false_split_upper_bracket=bound['bracket'],
                    exact_excess_p=str(excess),excess_rejected=bool(excess<=exact.ALPHA),
                    best_final100_upper=exact.upper(errors,100,exact.ALPHA)['value'],
                    final100_possible=bool(exact.cdf(errors,100,exact.TARGET)<=exact.ALPHA),
                    error_replicates=[r['replicate']for r in rows if r['null_error']]))
            assert all(a['seed']==b['seed']and a['case']==b['case']for a,b in
                zip(arms['production'],arms['pooled_prevalence']))
        primary=[r for r in summaries if r['arm']=='pooled_prevalence']
        assert len(primary)==2
        write_json(results,args.out/'cases.json')
        write_json(summaries,args.out/'performance.json')
        write_json(dict(complete_null_workflows=40,independent_base_seeds=10,
            candidate_rejected=any(r['excess_rejected']for r in primary),
            final100_null_requirements_possible=all(r['final100_possible']for r in primary),
            target_error_control_established=False,power_pending=True,
            full_pilot_complete=False,biological_specimens=0,
            same_prespecified_n10_look=True,additional_hypothesis_looks=0,
            scientific_criteria_changed=False,confidence_checks=checks),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU exact-rational planned null look; no inference rerun')
    except BaseException as error:
        finish_run(args.out,rec,error=error);raise


if __name__=='__main__':main()
