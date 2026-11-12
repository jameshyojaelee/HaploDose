"""Independent solver probe selected solely by recorded numerical instability."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logsumexp
from scipy.stats import norm

from haplodose.provenance import start_run, finish_run, file_record


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def independent_objective(parameters,z):
    means=np.asarray(parameters[:2]);sd=np.exp(parameters[2]);logit=parameters[3]
    weights=-np.logaddexp(0.,np.array([logit,-logit]))
    density=norm.logpdf(z,means,sd)+weights
    return -float(np.mean(logsumexp(density,axis=1)))


def canonical(parameters):
    p=np.array(parameters,copy=True)
    if p[0]>p[1]:p=p[[1,0,2,3]];p[3]*=-1
    return p


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    args=p.parse_args();args.out.mkdir(exist_ok=False)
    root=args.task.parents[1];prior=root/'dev/rna_population_baseline_20261109'
    baseline_path=prior/'source_v3/benchmarks/hardening/rna_population_baseline.py'
    generator_path=root/'dev/recovery_next_20261107/rank_source_v1/benchmarks/hardening/global_rank_generator.py'
    baseline=load('frozen_probe_baseline',baseline_path);generator=load('frozen_probe_generator',generator_path)
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',baseline=baseline_path,generator=generator_path)
    records=[]
    for part in ('calibration','evaluation'):
        directory=prior/(part+'_v2');manifest=json.loads((directory/'run_manifest.json').read_text())
        assert manifest['status']=='complete'
        record=manifest['output']['files']['records.jsonl'];assert file_record(record['path'])['sha256']==record['sha256']
        inputs[part+':manifest']=directory/'run_manifest.json';inputs[part+':records']=Path(record['path'])
        records.extend(json.loads(line) for line in Path(record['path']).read_text().splitlines())
    key=lambda row:(row['role'],row['replicate'])
    selected={key(r):r for r in sorted(records,key=lambda r:(-r['refinement_score_change'],*key(r)))[:4]}
    selected.update({key(r):r for r in sorted(records,key=lambda r:
        (-min(r['refinement_MAP_changes'],600-r['refinement_MAP_changes']),*key(r)))[:4]})
    for k,row in selected.items():
        assert file_record(row['arrays']['path'])['sha256']==row['arrays']['sha256']
        inputs[f'{k}:arrays']=Path(row['arrays']['path'])
    run=start_run(args.out,inputs,dict(selection='Top4 score changes plus top4 label-invariant MAP counts, role/replicate ties',
        solvers=['L-BFGS-B continuation','independent Nelder-Mead'],maxiter=10000),
        dict(scope='Numerical mechanism probe; no outcome-based model selection or scientific gate change'),mode='rna_population_optimizer_probe')
    try:
        results=[];summary=[]
        for k,row in sorted(selected.items()):
            condition,role,n=baseline.ROLES[row['role']]
            y,reference,_,_=generator.generate(np.random.default_rng(np.random.SeedSequence(row['seed_key'])),condition,600,240)
            assert generator.array_hash(y)==row['input_sha256']['counts']
            y=y.astype(float);order=np.lexsort(y[:,::-1].T);ys=y[order]
            x=np.log((ys+.5)/(ys.sum(axis=1)[:,None]*reference[None]+.5));x-=x.mean(axis=0)
            with np.load(row['arrays']['path']) as a:loading=a['loading'];old_map=a['probability'].argmax(axis=1)
            score=x@loading;variance=float(np.mean((score-score.mean())**2));assert variance==row['projection_variance']
            z=((score-score.mean())/np.sqrt(variance))[:,None];ll1=float(np.sum(-.5*(np.log(2*np.pi)+z[:,0]**2)))
            minimum=min(sum(float((value-np.mean(group))**2) for group in (v[:i],v[i:]) for value in group)
                for v in [np.sort(z[:,0])] for i in range(1,len(z)))/len(z)
            bounds=[(float(z.min()),float(z.max()))]*2+[(.5*np.log(minimum),0.),(None,None)]
            starting=[]
            for fit in row['fits']:
                starting.append(np.r_[fit['means'],np.log(fit['sd']),np.log(fit['weights'][1]/fit['weights'][0])])
            candidates=[];saved_probability=[]
            for label,initial,old in zip(('coarse','fine'),starting,row['fits']):
                value=independent_objective(initial,z)
                discrepancy=abs(max(0.,2*(-len(z)*value-ll1))-old['score'])
                budget=256*len(z)*np.finfo(float).eps*(1+abs(value)+abs(ll1)/len(z))
                assert discrepancy<=budget,(k,label,discrepancy,budget)
                parameters=canonical(initial);probability=np.exp(norm.logpdf(z,parameters[:2],np.exp(parameters[2]))
                    -np.logaddexp(0.,np.array([parameters[3],-parameters[3]])))
                probability/=probability.sum(axis=1,keepdims=True)
                restored=np.empty_like(probability);restored[order]=probability
                saved_probability.append(restored)
                if label=='fine':assert np.array_equal(restored.argmax(axis=1),old_map)
                for solver in ('L-BFGS-B','Nelder-Mead'):
                    if solver=='L-BFGS-B':
                        fitted=minimize(lambda p:baseline.mixture_objective(p,z),initial,method=solver,jac=True,bounds=bounds,
                            options=dict(maxiter=10000,ftol=np.finfo(float).eps,gtol=np.sqrt(np.finfo(float).eps),maxls=100))
                    else:
                        fitted=minimize(lambda p:independent_objective(p,z),initial,method=solver,bounds=bounds,
                            options=dict(maxiter=10000,xatol=np.sqrt(np.finfo(float).eps),fatol=np.finfo(float).eps))
                    target=independent_objective(fitted.x,z)
                    gradient=baseline.mixture_objective(fitted.x,z)[1]
                    actual=2*(-len(z)*target-ll1)
                    fitted_density=norm.logpdf(z,fitted.x[:2],np.exp(fitted.x[2]))-np.logaddexp(0.,np.array([fitted.x[3],-fitted.x[3]]))
                    fitted_labels=fitted_density.argmax(axis=1)
                    minimum_group=int(np.bincount(fitted_labels,minlength=2).min())
                    record=dict(role=k[0],replicate=k[1],start=label,solver=solver,success=bool(fitted.success),
                        message=str(fitted.message),iterations=int(fitted.nit),parameters=canonical(fitted.x).tolist(),
                        score=actual,score_gain_over_start=actual-old['score'],max_gradient=float(np.max(np.abs(gradient))),
                        independent_objective_difference=abs(target-fitted.fun),minimum_MAP_group=minimum_group,
                        material_supported_gain_over_saved_fine=bool(actual-row['score']>budget and minimum_group>=30))
                    results.append(record)
                    if fitted.success:candidates.append(record)
            assert candidates
            best=max(candidates,key=lambda r:r['score'])
            coarse,fine=saved_probability
            probability_change=min(float(np.mean(np.abs(coarse-fine))),float(np.mean(np.abs(coarse-fine[:,::-1]))))
            literal=row['refinement_MAP_changes'];aligned=min(literal,600-literal)
            coarse_map=coarse.argmax(axis=1);fine_map=fine.argmax(axis=1)
            if literal>300:coarse_map=1-coarse_map
            changed=coarse_map!=fine_map
            assert int(changed.sum())==aligned
            summary.append(dict(role=k[0],replicate=k[1],saved_coarse_score=row['coarse_score'],saved_fine_score=row['score'],
                fine_lost_score_vs_coarse=row['coarse_score']-row['score'],saved_literal_MAP_changes=literal,
                label_invariant_MAP_changes=aligned,aligned_mean_absolute_probability_change=probability_change,
                changed_cells_fine_median_max_probability=float(np.median(fine[changed].max(axis=1))) if changed.any() else None,
                saved_fine_minimum_MAP_group=row['minimum_MAP_group'],
                best_score=best['score'],best_gain_over_saved_fine=best['score']-row['score'],
                best_solver=best['solver'],best_start=best['start']))
        pd.DataFrame(results).to_csv(args.out/'solver_runs.tsv',sep='\t',index=False)
        pd.DataFrame(summary).to_csv(args.out/'cases.tsv',sep='\t',index=False)
        (args.out/'report.json').write_text(json.dumps(dict(selected_cases=len(selected),solver_runs=len(results),
            cases=summary,scientific_gates_changed=False,method_promoted=False,global_optimality_certified=False),indent=2)+'\n')
        for record in run['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(args.out,run,backend='cpu')
    except Exception as error:finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
