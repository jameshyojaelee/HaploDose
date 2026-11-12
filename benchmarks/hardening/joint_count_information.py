"""Postfit information diagnostic, not independent confirmation or calibration.

Compare exact shared-membership state inference on the full observed count with
its two halves. All cells start with the same proposal-fitted population mixing
probability. Never reuse proposal-derived cell membership as a full-count prior.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import expit
from haplodose.clones import CloneParams, _thin_diff_params
from haplodose.profile import DIFF_LEVELS, DIFF_PRIOR, credible_difference_runs
from haplodose.provenance import start_run, finish_run, file_record


def load(name):
    path=Path(__file__).with_name(name+'.py')
    spec=importlib.util.spec_from_file_location('information_'+name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True);p.add_argument('--previous-oracle',type=Path,required=True)
    p.add_argument('--previous-probe',type=Path)
    a=p.parse_args();a.out.mkdir(exist_ok=False);model=load('shared_membership_hmm')
    fixture=load('membership_centering_factorial');root=Path(__file__).parent
    old=json.loads((a.previous_oracle/'run_manifest.json').read_text())
    if old['status']!='complete':raise ValueError('Original diagnostic incomplete')
    inputs={k:Path(r['path']) for k,r in old['input']['files'].items()}
    inputs.update(script=Path(__file__),model=root/'shared_membership_hmm.py',
                  independent_reference=root/'membership_centering_factorial.py',
                  previous_manifest=a.previous_oracle/'run_manifest.json',
                  previous_result=a.previous_oracle/'summary.json')
    if a.previous_probe is not None:
        inputs['previous_probe']=a.previous_probe/'summary.json'
    for r in old['input']['files'].values():
        if file_record(r['path'])['sha256']!=r['sha256']:raise ValueError('Frozen source changed')
    rec=start_run(a.out,inputs,dict(seeds=[20261108,20261109,20261110],mixing_weight='proposal-fitted population proportion'),
        dict(scope='Selected-event conditional information diagnostic; no independent confirmation, no error or recovery claim',
             estimand='relative RNA-state trajectory and latent memberships conditional on fixed nuisance/model'),
        mode='joint_count_information')
    try:
        (a.out/'reference_checks.json').write_text(json.dumps(fixture.checks(model),indent=2))
        previous=json.loads(inputs['previous_result'].read_text());rows=[];identities=[]
        previous_probe=None if a.previous_probe is None else json.loads(inputs['previous_probe'].read_text())['cases']
        pp=_thin_diff_params(CloneParams(split_mode='thin').profile)
        for seed in (20261108,20261109,20261110):
            events=[x for x in json.loads(inputs[f'{seed}:trace.json'].read_text()) if x['stage']=='event_fit']
            if len(events)!=1:raise ValueError('Locked event population changed')
            with np.load(inputs[f'{seed}:event{events[0]["index"]}']) as z:
                O,T,E,r,q,cols,parameters=[z[k] for k in ('O','Ytest','E','r','q','cols','parameters')]
            chrom=pd.read_csv(inputs[f'{seed}:fit/bins.tsv'],sep='\t').chromosome.to_numpy()
            beta=float(pd.read_csv(inputs[f'{seed}:fit/specimens.tsv'],sep='\t').set_index('specimen').loc['S1','beta'])
            population=float(expit(parameters[2]));prior=np.full(len(O),population)
            boundary=dict(left_neutral=bool(cols[0]>0 and chrom[cols[0]-1]==chrom[cols[0]]),
                          right_neutral=bool(cols[-1]+1<len(chrom) and chrom[cols[-1]+1]==chrom[cols[-1]]))
            def build(y,e,size,probability):
                return model.build(y,e,size,probability,DIFF_LEVELS*beta*np.log(1.5),DIFF_PRIOR,pp.p_change,**boundary)
            models=dict(proposal=build(O,E,r,prior),test=build(T,E,r,prior),
                        joint=build(O+T,2*E,2*r,prior),frozen_conditional=build(T,E,r,q))
            # Under this NB working law, fold scaling preserves E/r. The
            # null-relative component log likelihoods add before integrating
            # the one shared membership and state trajectory. Priors enter once.
            sequential=dict(models['proposal'],ll=models['proposal']['ll']+models['test']['ll'])
            eps=np.finfo(float).eps;nops=64+4*O.size
            gamma=nops*eps/(1-nops*eps)
            allowance=gamma*max(1.,float(np.abs(models['joint']['ll']).sum()))
            discrepancy=float(np.max(np.abs(models['joint']['ll']-sequential['ll'])))
            if discrepancy>allowance:raise ArithmeticError('NB fold composition identity failed')
            fits={name:model.enumerate_states(value) for name,value in models.items()}
            sequential_fit=model.enumerate_states(sequential)
            posterior_error=float(np.max(np.abs(fits['joint']['posterior']-sequential_fit['posterior'])))
            if posterior_error>4*allowance:raise ArithmeticError('Sequential/joint marginalization disagree')
            if abs(fits['joint']['log_evidence']-sequential_fit['log_evidence'])>allowance:
                raise ArithmeticError('Joint evidence differs from sequential composition')
            original=next(x for x in previous if x['seed']==seed and x['membership']=='soft')
            np.testing.assert_array_equal(fits['frozen_conditional']['posterior'],original['posterior'])
            if fits['frozen_conditional']['log_evidence']!=original['log_evidence']:
                raise AssertionError('Original conditional result changed')
            identities.append(dict(seed=seed,component_error=discrepancy,posterior_error=posterior_error,
                                   accumulation_allowance=allowance,original_conditional_exact=True))
            for name,fit in fits.items():
                reset=np.r_[True,np.zeros(len(cols)-1,bool)]
                cores=[[int(cols[x]),int(cols[y]),int(s)] for x,y,s in
                       credible_difference_runs(fit['path'],fit['posterior'],reset,pp)]
                # A numerical ambiguity never changes the scientific criterion.
                support=fit['posterior'][np.arange(len(cols)),fit['path']]
                margin=float(np.min(np.abs(support-pp.posterior)))
                if previous_probe is not None:
                    old_fit=next(x for x in previous_probe if x['seed']==seed and x['method']==name)
                    np.testing.assert_array_equal(fit['posterior'],old_fit['posterior'])
                    if cores!=old_fit['cores'] or fit['log_evidence']!=old_fit['log_evidence']:
                        raise AssertionError('Reporting correction changed scientific results')
                row=dict(seed=seed,method=name,population_proportion=population,log_evidence=fit['log_evidence'],
                    posterior=fit['posterior'].tolist(),cores=cores,posterior_rule=pp.posterior,min_bins=pp.min_bins,
                    path=fit['path'].tolist(),minimum_decoded_state_posterior=float(support.min()),
                    numerical_boundary_metric='decoded directional-state posterior used by credible_difference_runs',
                    minimum_non_neutral_posterior=float(np.min(1-fit['posterior'][:,1])),
                    numerical_boundary_unresolved=bool(margin<=4*allowance),
                    interpretation='Adaptively selected interval statistic, not an independently confirmed CNA event')
                rows.append(row);print(json.dumps({k:v for k,v in row.items() if k!='posterior'}),flush=True)
        if len(rows)!=12:raise AssertionError('Missing locked comparisons')
        (a.out/'summary.json').write_text(json.dumps(dict(cases=rows,identities=identities,
            decision='No promotion or global error claim; any useful joint procedure requires full-workflow null/power tests'),indent=2))
        baseline=json.loads((a.task/'baseline/source_manifest.json').read_text());core=0
        for relative,digest in baseline.items():
            if file_record(a.task/'baseline'/relative)['sha256']!=digest:raise AssertionError('Baseline changed')
            if relative.startswith(('python/','tests/','crates/')):
                if file_record(a.task.parents[1]/relative)['sha256']!=digest:raise AssertionError('Core changed')
                core+=1
        for r in rec['input']['files'].values():
            if file_record(r['path'])['sha256']!=r['sha256']:raise AssertionError('Frozen input changed')
        (a.out/'preservation.json').write_text(json.dumps(dict(baseline_files=len(baseline),core_files=core,
            input_hashes=len(rec['input']['files']),original_conditional_reproductions=3),indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
