"""Oracle-only expected-score attribution after frozen RNA prediction."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from haplodose.provenance import start_run,finish_run,file_record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);panel=a.task/'panel_v1'
    helper_path=Path(__file__).with_name('reference_rate_prediction.py')
    spec=importlib.util.spec_from_file_location('expected_reference_objective',helper_path)
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    means=a.task.parent/'profile_pool_mechanism_20261109/panel_v1/primary_means_and_profiles.npz'
    inputs=dict(script=Path(__file__),helper=helper_path,design=a.task/'TARGET_DESIGN.md',
        predictions=panel/'predictions.npz',frozen=panel/'prediction_freeze.json',means=means,
        audit=a.task/'audit_v1/run_manifest.json')
    assert json.loads(inputs['audit'].read_text())['status']=='complete'
    for record in json.loads(inputs['frozen'].read_text()).values():assert file_record(record['path'])['sha256']==record['sha256']
    run=start_run(a.out,inputs,{},dict(scope='Post hoc oracle expected-gradient attribution; not a candidate fit'),
        mode='reference_rate_target')
    try:
        with np.load(inputs['predictions']) as z:v={key:z[key] for key in z.files}
        with np.load(means) as z:M=z['used_gene_mean']
        rows=[];stored={}
        for group in ('type0','type1'):
            true=M[v[group+'__fit']];depth=true.sum(axis=1)
            pooled=true.sum(axis=0)/depth.sum()
            equal_depth=np.full(len(depth),depth.mean());equal_means=equal_depth[:,None]*(true/depth[:,None])
            conditions=dict(generating=(true,depth),common_rate=(depth[:,None]*pooled,depth),
                            equal_depth=(equal_means,equal_depth))
            r=1/v[group+'__phi']
            for condition,(Y,L) in conditions.items():
                q=Y.sum(axis=0)/L.sum();z=np.log(q[:-1]/q[-1])
                _,gradient=helper.objective(z,Y,L,r,q)
                g=np.array([math.fsum(float(r[j])*(float(d)*float(q[j])-float(y))/(float(r[j])+float(d)*float(q[j]))
                                     for d,y in zip(L,Y[:,j])) for j in range(len(q))])
                independent=(g-q*math.fsum(g))/math.fsum(Y.ravel());independent=independent[:-1]
                complex_gradient=[]
                for j in range(len(z)):
                    point=z.astype(complex);point[j]+=1e-25j
                    complex_gradient.append(helper.scalar_complex_objective(point,Y,L,r,q).imag/1e-25)
                operations=64*Y.size;u=np.finfo(float).eps
                budget=operations*u/(1-operations*u)*max(1,float(np.max(np.abs(gradient))))
                discrepancy=max(float(np.max(np.abs(gradient-independent))),
                                float(np.max(np.abs(gradient-np.array(complex_gradient)))) )
                assert discrepancy<=budget
                maximum=float(np.max(np.abs(gradient)))
                if condition!='generating':assert maximum<=budget
                # Poisson score limit requires no chosen small dispersion.
                poisson=np.array([math.fsum(float(d)*float(q[j])-float(y) for d,y in zip(L,Y[:,j]))
                                  for j in range(len(q))])
                poisson=(poisson-q*math.fsum(poisson))/math.fsum(Y.ravel())
                assert np.max(np.abs(poisson))<=budget
                rows.append(dict(group=group,condition=condition,max_projected_expected_NB_gradient=maximum,
                    independent_gradient_discrepancy=discrepancy,operation_budget=budget,
                    resolved_nonzero=bool(maximum>budget),max_expected_Poisson_gradient=float(np.max(np.abs(poisson)))))
                stored[group+'__'+condition+'__gradient']=gradient
        np.savez_compressed(a.out/'expected_gradients.npz',**stored)
        (a.out/'report.json').write_text(json.dumps(dict(conditions=rows,oracle_only=True,
            post_hoc=True,no_predictions_or_selection_changed=True,production_changed=False,
            scientific_criteria_changed=False,biological_specimens=0),indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out,run,error=error);raise
    finish_run(a.out,run,backend='CPU independent analytic/scalar/complex expected scores')


if __name__=='__main__':main()
