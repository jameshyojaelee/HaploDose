"""Independent high-precision forward rates and calibration/contrast recount."""
import argparse
from decimal import Decimal,localcontext
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record


def scalar_rates(depth,gamma,effect=None):
    with localcontext() as context:
        context.prec=80
        ratio=Decimal.from_float(float(depth))/Decimal(9000);root=ratio.sqrt();weights=[]
        for j,g in enumerate(gamma):
            value=1/root if g==-.5 else root if g==.5 else Decimal(1)
            if effect is not None:value*=Decimal.from_float(float(effect[j])).exp()
            weights.append(value)
        total=sum(weights,Decimal(0))
        return [weight/total for weight in weights]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);panel=a.task/'panel_v1'
    inputs=dict(script=Path(__file__),manifest=panel/'run_manifest.json',models=panel/'forward_models.npz',
        contrasts=panel/'contrast_errors.tsv',report=panel/'report.json')
    previous=json.loads(inputs['manifest'].read_text());assert previous['status']=='complete'
    for group in ('input','output'):
        for record in previous[group]['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
    run=start_run(a.out,inputs,{},dict(scope='Independent deterministic forward-model audit, no biological inference'),
        mode='depth_reference_anchor_audit')
    try:
        with np.load(inputs['models']) as z:v={key:z[key] for key in z.files}
        frame=pd.read_csv(inputs['contrasts'],sep='\t',float_precision='round_trip',keep_default_na=False)
        assert set(frame.event)=={'null','gain'}
        assert len(frame)==24 and not frame[['condition','event','method']].duplicated().any()
        delta=.6*math.log(1.5);rows=[];mean_checks=0;gradient_checks=[]
        for condition in ('matched','shifted','reference_spread','zero_slope'):
            reference=v[condition+'__reference_depth'];query=v[condition+'__query_depth'];gamma=v[condition+'__gamma']
            mu=v[condition+'__reference_rate'];b=v[condition+'__b'];origin=float(v[condition+'__origin'])
            u=np.finfo(float).eps;operations=128*len(mu)*(len(query)+len(reference));budget=operations*u/(1-operations*u)
            exact_reference=np.zeros(40)
            for depth in reference:exact_reference+=np.array([float(q) for q in scalar_rates(depth,gamma)])*float(depth)
            exact_reference/=math.fsum(reference)
            assert (np.abs(exact_reference-mu)<=budget*np.maximum(1,np.abs(mu))).all()
            # Scalar calibration equations and Hessian independently reconstruct
            # pooled reference rates and their local convexity.
            predicted=np.zeros(40);H=np.zeros((40,40))
            for depth in reference:
                weights=[math.exp(float(bb)+float(g)*(math.log(float(depth))-origin)) for bb,g in zip(b,gamma)]
                q=np.array([x/math.fsum(weights) for x in weights]);w=float(depth)/math.fsum(reference)
                predicted+=w*q;H+=w*(np.diag(q)-np.outer(q,q))
            residual=float(np.max(np.abs(predicted-mu)));smallest=float(np.linalg.eigvalsh(H[:-1,:-1]).min())
            assert smallest>0
            solvers=json.loads(inputs['report'].read_text())['solver_records']
            record=[r for r in solvers if r['condition']==condition and r['tolerance']==1e-12 and not r.get('origin_shift')][0]
            assert residual<=max(1e-12,record['operation_budget'])+budget
            gradient_checks.append(dict(condition=condition,reference_residual=residual,minimum_anchored_Hessian_eigenvalue=smallest))
            for event in ('null','gain'):
                effect=np.r_[np.full(10,delta),np.zeros(30)] if event=='gain' else np.zeros(40)
                exact=np.zeros(40)
                for depth in query:exact+=np.array([float(q) for q in scalar_rates(depth,gamma,effect)])*float(depth)
                observed=v[condition+'__'+event+'__mean']
                assert (np.abs(observed-exact)<=budget*np.maximum(1,np.abs(exact))).all();mean_checks+=len(exact)
                for method in ('current','origin_shift_only','anchored'):
                    E=v[condition+'__'+method+'__expected']
                    if method=='current':assert (np.abs(E-query.sum()*mu)<=budget*np.maximum(1,np.abs(E))).all()
                    logratio=np.array([math.log(float(o)/float(e)) for o,e in zip(observed,E)])
                    ordered=sorted(logratio);median=(ordered[19]+ordered[20])/2
                    residuals=logratio-median
                    metrics=dict(max_abs_centred_log_error=max(abs(x-y) for x,y in zip(residuals,effect)),
                        first_region_mean=math.fsum(residuals[:10])/10,last_region_mean=math.fsum(residuals[-10:])/10)
                    row=frame.loc[frame.condition.eq(condition)&frame.event.eq(event)&frame.method.eq(method)].iloc[0]
                    for name,value in metrics.items():assert abs(row[name]-value)<=budget*(1+abs(value))
                    rows.append(dict(condition=condition,event=event,method=method,**metrics))
        # Independent analytic shifted-depth contrast including true gain sign.
        shifted=frame.loc[frame.condition.eq('shifted')&frame.method.eq('current')]
        for row in shifted.itertuples():
            target=-math.log(2)+(delta if row.event=='gain' else 0)
            assert abs(row.first_region_mean-target)<=128*np.finfo(float).eps*(1+abs(target))
        (a.out/'report.json').write_text(json.dumps(dict(contrast_rows_recounted=len(rows),
            decimal80_query_bin_means_checked=mean_checks,calibration_checks=gradient_checks,
            shifted_depth_analytic_contrast_verified=True,scientific_criteria_unchanged=True,
            biological_validation=False),indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out,run,error=error);raise
    finish_run(a.out,run,backend='CPU Decimal80/scalar calibration and contrast references')


if __name__=='__main__':main()
