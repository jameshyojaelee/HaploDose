"""Conditional calibration-bank uncertainty, separate from the frozen holdout gate."""
import argparse
import importlib.util
import json
from pathlib import Path
from haplodose.provenance import start_run,finish_run,file_record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    source=a.task/'summary_source_v1/benchmarks/hardening/joint_fresh_summary.py'
    spec=importlib.util.spec_from_file_location('frozen_bounds',source)
    stats=importlib.util.module_from_spec(spec);spec.loader.exec_module(stats)
    inputs=dict(script=Path(__file__),numerics=source,summary=a.task/'summary_v1/summary.json',
        audit=a.task/'audit_v1/report.json',audit_manifest=a.task/'audit_v1/run_manifest.json',
        design=a.task/'DESIGN.md')
    audit=json.loads(inputs['audit_manifest'].read_text())
    if audit['status']!='complete' or file_record(inputs['audit'])['sha256']!=audit['output']['files']['report.json']['sha256']:
        raise AssertionError('Fresh audit changed or incomplete')
    rec=start_run(a.out,inputs,dict(calibration_draws=199,reference_family=153),
        dict(scope='Distributional rank calculation only; empirical gate failures retained'),mode='rank_bank_precision')
    try:
        if stats.FAMILY!=153:raise AssertionError('Changed numerical confidence allocation')
        rows=[]
        for alpha,rank in ((.05,10),(.025,5),(.01,2),(.005,1)):
            # With n continuous exchangeable calibration scores, the probability
            # beyond the kth-largest score has Beta(k,n+1-k) distribution.
            # Its upper quantile solves P[Binomial(n,theta)<=k-1]=delta.
            upper=stats.bound(rank-1,199,False)
            lower,high=stats.independent_enclosure(rank-1,199,False)
            if not lower<=upper<=high:raise AssertionError('Independent bank bound differs')
            rows.append(dict(rank_gate=alpha,largest_score_rank=rank,
                marginal_rank_error_bound=rank/200,conditional_tail_probability_upper=upper,
                numerical_enclosure=[lower,high],tail_probability=.05/153,
                assumption='Calibration and future scores exchangeable under exactly the same complete workflow/null law; ties conservative'))
        (a.out/'report.json').write_text(json.dumps(dict(rows=rows,
            empirical_failure_overridden=False,new_inference_draws=0,new_biological_units=0,
            interpretation='An unconditional rank bound is not a fixed-bank conditional bound. These '
                'model-conditional calculations do not replace the prespecified independent holdout criterion '
                'or supply biological validity.'),indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
