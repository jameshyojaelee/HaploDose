"""High-precision heldout score and explicit information-leakage replay."""
import argparse
from decimal import Decimal,localcontext
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record


COEFFICIENTS={}


def decimal_score(Y,L,q,r,digits):
    with localcontext() as context:
        context.prec=digits
        terms=[]
        for i in range(len(Y)):
            for j in range(Y.shape[1]):
                y=int(Y[i,j]);assert y==Y[i,j]
                shape=Decimal.from_float(float(r[j]));key=(digits,float(r[j]))
                # For integer y the gamma ratio is a rising-factorial product:
                # Gamma(r+y)/(Gamma(r)*y!)=product_k (r+k)/(k+1).
                # This avoids both SciPy's gammaln and any extra dependency.
                products=COEFFICIENTS.setdefault(key,[Decimal(1)])
                while len(products)<=y:
                    k=len(products)-1;products.append(products[-1]*(shape+k)/Decimal(k+1))
                mean=Decimal.from_float(float(L[i]))*Decimal.from_float(float(q[j]))
                logtotal=(shape+mean).ln()
                terms.append(products[y].ln()+shape*(shape.ln()-logtotal)+y*(mean.ln()-logtotal))
        return sum(terms,Decimal(0))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);panel=a.task/'panel_v1'
    source=Path(__file__).with_name('reference_rate_prediction.py')
    inputs=dict(script=Path(__file__),helper=source,manifest=panel/'run_manifest.json',
        predictions=panel/'predictions.npz',freeze=panel/'prediction_freeze.json',
        scores=panel/'heldout_scores.tsv',truth_errors=panel/'post_freeze_truth_errors.tsv',
        selection=panel/'selection_before_truth.json',optimizations=panel/'optimizations.json',checks=panel/'numerical_checks.json')
    previous=json.loads(inputs['manifest'].read_text());assert previous['status']=='complete'
    for group in ('input','output'):
        for record in previous[group]['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
    for record in json.loads(inputs['freeze'].read_text()).values():assert file_record(record['path'])['sha256']==record['sha256']
    run=start_run(a.out,inputs,{},dict(scope='High-precision score and information-isolation validation, no biological claim'),
        mode='reference_rate_prediction_audit')
    try:
        with np.load(inputs['predictions']) as z:v={key:z[key] for key in z.files}
        prior=a.task.parent/'profile_pool_export_20261109/parity_v1/fitted_new'
        Y=np.load(prior/'bin_counts.npy').astype(float);meta=pd.read_csv(prior/'cells.tsv.gz',sep='\t')
        frame=pd.read_csv(inputs['scores'],sep='\t',float_precision='round_trip')
        scores=[];numerical=[]
        for group in ('type0','type1'):
            fit=meta.role.eq('reference')&meta.ref_group.eq(group)&meta.ref_split.eq('fit')
            hold=meta.role.eq('reference')&meta.ref_group.eq(group)&meta.ref_split.eq('holdout')
            assert np.array_equal(fit,v[group+'__fit']) and np.array_equal(hold,v[group+'__holdout'])
            fitY=Y[fit];L=fitY.sum(axis=1)
            total=np.array([math.fsum(fitY[:,j]) for j in range(Y.shape[1])]);total=np.maximum(total,.1)
            pooled=total/math.fsum(total)
            equal=np.array([math.fsum(float(fitY[i,j])/float(L[i]) for i in range(len(fitY)))/len(fitY)
                            for j in range(Y.shape[1])]);equal/=math.fsum(equal)
            u=np.finfo(float).eps;operations=32*fitY.size;gamma=operations*u/(1-operations*u)
            for label,reference in [('pooled',pooled),('equal',equal)]:
                assert (np.abs(reference-v[group+'__'+label])<=gamma*np.maximum(1,np.abs(reference))).all()
            phi=[]
            for j,q in enumerate(pooled):
                means=[float(l)*q for l in L]
                numerator=math.fsum((float(y)-m)**2-m for y,m in zip(fitY[:,j],means))
                denominator=math.fsum(m*m for m in means)
                phi.append(min(10.,max(.001,numerator/denominator)))
            phi=np.array(phi);assert (np.abs(phi-v[group+'__phi'])<=gamma*np.maximum(1,np.abs(phi))).all()
            r=1/v[group+'__phi'];held=Y[hold];LH=held.sum(axis=1)
            for method in ('pooled','equal','nb','nb_1e-08'):
                q=v[group+'__'+method];assert (q>0).all() and abs(math.fsum(q)-1)<=gamma
                lo=decimal_score(held,LH,q,r,60);hi=decimal_score(held,LH,q,r,80)
                difference=float(abs(lo-hi));assert difference<10.**-50
                score=float(hi)
                if method!='nb_1e-08':
                    recorded=float(frame.loc[frame.group.eq(group)&frame.method.eq(method),'holdout_composite_log_score'].iloc[0])
                    discrepancy=abs(score-recorded)
                    assert discrepancy<=gamma*max(1,abs(score))
                    numerical.append(dict(group=group,method=method,decimal60_decimal80_discrepancy=difference,
                        independent_score_discrepancy=discrepancy,operation_budget=gamma*max(1,abs(score))))
                scores.append(dict(group=group,method=method,decimal80_score=score))
        pd.DataFrame(scores).to_csv(a.out/'high_precision_scores.tsv',sep='\t',index=False)
        selections=json.loads(inputs['selection'].read_text())
        for selection in selections:
            f=[r for r in scores if r['group']==selection['group']]
            fine=sorted([r for r in f if r['method']!='nb_1e-08'],key=lambda r:-r['decimal80_score'])
            assert [r['method'] for r in fine]==selection['fine_ranking']
            assert selection['refinement_choice_stable']
            budget=max(r['operation_budget'] for r in numerical if r['group']==selection['group'])
            assert fine[0]['decimal80_score']-fine[1]['decimal80_score']>2*budget
        original=Path(previous['input']['files']['means']['path'])
        with np.load(original) as z:M=z['used_gene_mean']
        # Alter the latent truth in a fresh file, never a frozen measurement.
        perturb=M*np.exp(np.linspace(-.5,.5,M.shape[1]))[None]
        with (a.out/'perturbed_truth.npz').open('xb') as stream:np.savez_compressed(stream,used_gene_mean=perturb)
        errors=pd.read_csv(inputs['truth_errors'],sep='\t',float_precision='round_trip');checked=0
        for row in errors.itertuples():
            mask=v[row.group+'__'+row.target_split]
            target=np.array([math.fsum(M[mask,j]) for j in range(M.shape[1])]);target/=math.fsum(target)
            q=v[row.group+'__'+row.method];x=np.log(q);x-=np.median(x);y=np.log(target);y-=np.median(y)
            actual=math.fsum((a-b)**2 for a,b in zip(x,y))/len(x)
            assert abs(actual-row.centred_log_rate_MSE)<=gamma*(1+abs(actual));checked+=1
        for label,means in [('original',original),('perturbed',a.out/'perturbed_truth.npz')]:
            destination=a.out/f'replay_{label}'
            result=subprocess.run([sys.executable,'-I',str(source),'--task',str(a.task),'--out',str(destination),
                '--means',str(means)],capture_output=True,text=True,timeout=1800,cwd=a.out)
            (a.out/f'{label}.stdout').write_text(result.stdout);(a.out/f'{label}.stderr').write_text(result.stderr)
            if result.returncode:raise RuntimeError(f'{label} replay failed with {result.returncode}')
            with np.load(destination/'predictions.npz') as z:
                assert set(z.files)==set(v)
                for key in z.files:assert np.array_equal(z[key],v[key]),key
            for name in ('heldout_scores.tsv','selection_before_truth.json','optimizations.json','numerical_checks.json'):
                assert (destination/name).read_bytes()==(panel/name).read_bytes(),name
        assert (a.out/'replay_original/post_freeze_truth_errors.tsv').read_bytes()==inputs['truth_errors'].read_bytes()
        assert (a.out/'replay_perturbed/post_freeze_truth_errors.tsv').read_bytes()!=inputs['truth_errors'].read_bytes()
        (a.out/'report.json').write_text(json.dumps(dict(reference_types=2,heldout_score_rows=8,
            high_precision_digits=[60,80],score_checks=numerical,truth_error_checks=checked,
            full_original_and_perturbed_replays_exact=True,predictions_and_selection_isolated_from_truth=True,
            scientific_criteria_changed=False,biological_validation=False),indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out,run,error=error);raise
    finish_run(a.out,run,backend='CPU independent scalar/Decimal and complete analysis replays')


if __name__=='__main__':main()
