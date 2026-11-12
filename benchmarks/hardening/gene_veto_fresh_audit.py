"""Independent fresh-workflow event recount, generator correspondence and bounds."""
import argparse
import itertools
import json
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import brentq
from scipy.special import logsumexp
from scipy.stats import beta, binom
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record
from haplodose.simulate import simulate


CONDITIONS = ('matched','state','two_clone')
MASTER = 20261109


def read_counts(path):
    with h5py.File(path) as h:
        ids=h['cell_ids'].asstr()[:];genes=h['genes'].asstr()[:]
        X=sparse.csr_matrix((h['X/data'][:],h['X/indices'][:],h['X/indptr'][:]),shape=(len(ids),len(genes)))
    return X,ids,genes


def source_reference(rep,phase):
    captured=[]
    def observe(frame,event,value):
        if event=='return' and frame.f_code.co_name=='draw' and frame.f_code.co_filename==simulate.__code__.co_filename:
            loc=frame.f_locals
            if loc['cols']['specimen']=='S1':captured.append((loc['L'].copy(),loc['share'].copy()))
    assert sys.getprofile() is None
    sys.setprofile(observe)
    try:
        X,genes,meta,pos,events,_=simulate(seed=MASTER*10000+phase*1000+rep,
            specimens=[dict(name='S1',n_cells=600,normal_frac=0.,clones=[dict(frac=1.,events=[
                ('chr3',0.,1.,3),('chr8',.5,1.,1)])])])
    finally:sys.setprofile(None)
    assert len(captured)==1 and set(events.chromosome)=={'chr3','chr8'}
    assert set(events.cn)=={1,3}
    rows=np.flatnonzero(meta.specimen.eq('S1'));L,share=captured[0]
    return X[:rows[0]],L[:,None]*share,meta.cell_id.to_numpy(),np.array(genes)


def independent_score(folder, truth):
    cells=pd.read_csv(folder/'cells.tsv.gz',sep='\t')
    assert cells.cell_id.is_unique
    cells=cells[cells.specimen.eq('S1')].merge(truth[['cell_id','truth_clone']],on='cell_id',validate='one_to_one')
    assert len(cells)==600 and cells.truth_clone.notna().all()
    table=pd.read_csv(folder/'clones.tsv',sep='\t')
    names=sorted(table.loc[table.specimen.eq('S1') & ~table.normal_like,'clone'])
    diff=pd.read_csv(folder/'differences.tsv',sep='\t')
    chromosomes=sorted(set(diff.loc[diff.specimen.eq('S1') & diff.clone_a.isin(names) & diff.clone_b.isin(names),'chromosome']))
    populations=sorted(set(cells.truth_clone))
    n=min(len(populations),len(names));best=0
    # Exhaustive matching independent of the producer's Hungarian algorithm.
    for cols in itertools.permutations(names,n):
        for pops in itertools.permutations(populations,n):
            best=max(best,sum(int(sum((cells.truth_clone==pop)&(cells.clone==col)&cells.malignant_call))
                for pop,col in zip(pops,cols)))
    return dict(reported_groups=len(names),null_error=len(names)>1,
        exact_chr12_recovery=len(names)==2 and chromosomes==['chr12'],any_two_groups=len(names)>=2,
        difference_chromosomes=chromosomes,assignment_correct=best,assignment_n=600,tumour_calls=int(sum(cells.malignant_call)))


def bound(k,n,upper,alpha):
    if upper:
        b=1. if k==n else float(beta.ppf(1-alpha,k+1,n-k))
        q=1. if k==n else brentq(lambda p:binom.cdf(k,n,p)-alpha,0.,1.,xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps)
    else:
        b=0. if k==0 else float(beta.ppf(alpha,k,n-k+1))
        q=0. if k==0 else brentq(lambda p:binom.sf(k-1,n,p)-alpha,0.,1.,xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps)
    assert abs(b-q)<=128*np.finfo(float).eps
    return b,q


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--pilot',action='store_true')
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    directories=[a.task/'pilot_v1'] if a.pilot else sorted(a.task.glob('chunk_*_v1'))
    assert len(directories)==(1 if a.pilot else 12)
    inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md')
    records=[]
    for folder in directories:
        manifest=json.loads((folder/'run_manifest.json').read_text());assert manifest['status']=='complete'
        inputs[str(folder.relative_to(a.task))+':manifest']=folder/'run_manifest.json'
        inputs[str(folder.relative_to(a.task))+':results']=folder/'results.json'
        records.extend(json.loads((folder/'results.json').read_text()))
    for row in records:
        case=Path(row['case']);key=f"{row['condition']}:{row['rep']}:{row['arm']}"
        for file in ('counts.h5','truth_cells.tsv','model_cells.tsv','positions.tsv'):
            inputs[key+':input:'+file]=case/'input'/file
        inputs[key+':generation']=case/'generation.npz'
        for file in ('cells.tsv.gz','clones.tsv','differences.tsv','qc.json','run_manifest.json'):
            inputs[key+':fit:'+file]=case/row['arm']/'fit'/file
        inputs[key+':trace']=case/row['arm']/'trace.json'
    if a.pilot:
        original=a.task/'original_v1/replay'
        prior=a.task.parent/'gene_veto_aligned_20261109/candidate_v1/replay'
        for name,folder in (('original',original),('prior',prior)):
            inputs[name+':summary']=folder/'panel/summary.json';inputs[name+':trace']=folder/'trace.json'
            for item in json.loads((folder/'trace.json').read_text()):inputs[name+':fold:'+str(item['model'])]=Path(item['path'])
            for seed in (20261031,20261107,20261108,20261109,20261110):
                for file in ('cells.tsv.gz','differences.tsv'):
                    inputs[f'{name}:{seed}:{file}']=folder/f'panel/seed_{seed}/fit'/file
    rec=start_run(a.out,inputs,dict(pilot=a.pilot),dict(scope='Fresh conditional simulation, not biological replication'),mode='gene_veto_fresh_audit')
    try:
        if a.pilot:
            current=json.loads(inputs['original:summary'].read_text());previous=json.loads(inputs['prior:summary'].read_text())
            for b,c in zip(previous['cases'],current['cases']):
                assert b['seed']==c['seed'] and b['criteria']==c['criteria'] and b['original_assertion_passed']==c['original_assertion_passed']
                assert c['original_assertion_passed']
                for file in ('cells.tsv.gz','differences.tsv'):
                    pd.testing.assert_frame_equal(pd.read_csv(inputs[f"original:{c['seed']}:{file}"],sep='\t'),
                        pd.read_csv(inputs[f"prior:{c['seed']}:{file}"],sep='\t'),check_exact=True)
            for c,b in zip(json.loads(inputs['original:trace'].read_text()),json.loads(inputs['prior:trace'].read_text())):
                with np.load(c['path']) as x,np.load(b['path']) as y:
                    assert set(x.files)==set(y.files)
                    for key in x.files:assert np.array_equal(x[key],y[key])
        n=2 if a.pilot else 100;phase=0 if a.pilot else 1
        expected=set(itertools.product(CONDITIONS,range(n),('comparator','candidate')))
        observed={(r['condition'],r['rep'],r['arm']) for r in records}
        assert observed==expected and len(observed)==len(records)
        checked_cases=set();paired={};audited=[]
        for row in records:
            condition,rep,arm=row['condition'],row['rep'],row['arm'];case=Path(row['case'])
            assert row['phase']==phase and row['seed']==MASTER*10000+phase*1000+rep
            truth=pd.read_csv(case/'input/truth_cells.tsv',sep='\t')
            model=pd.read_csv(case/'input/model_cells.tsv',sep='\t')
            assert not any(c.startswith('truth_') for c in model)
            pd.testing.assert_frame_equal(model,truth.drop(columns=[c for c in truth if c.startswith('truth_')]),check_exact=True)
            key=(condition,rep)
            if key not in checked_cases:
                X,ids,genes=read_counts(case/'input/counts.h5')
                assert np.array_equal(ids,truth.cell_id) and np.all(X.data>=0) and np.issubdtype(X.dtype,np.integer)
                with np.load(case/'generation.npz') as z:
                    query=np.flatnonzero(truth.specimen.eq('S1'));assert np.array_equal(query,z['query_rows']) and len(query)==600
                    delta=.6*np.log(1.5);sigma=delta*np.sqrt(.6*.4)
                    assert float(z['delta'])==delta and float(z['sigma'])==sigma
                    state=np.random.default_rng(np.random.SeedSequence([MASTER,phase,rep,1])).normal(0.,sigma,600)
                    labels=np.zeros(600,int);labels[np.random.default_rng(np.random.SeedSequence([MASTER,phase,rep,2])).permutation(600)[:240]]=1
                    assert np.array_equal(state,z['state']) and np.array_equal(labels,z['latent_labels'])
                    t=labels if condition=='two_clone' else np.zeros(600,int)
                    assert np.array_equal(t,z['truth']) and np.array_equal(np.where(t==0,'clone1','clone2'),truth.iloc[query].truth_clone)
                    pos=pd.read_csv(case/'input/positions.tsv',sep='\t');assert np.array_equal(pos.gene,genes)
                    region=pos.chromosome.eq('chr12').to_numpy() & (pos.start.to_numpy()<=10000000)
                    assert np.array_equal(region,z['region'])
                    effect=np.zeros(600) if condition=='matched' else state.copy()
                    if condition=='two_clone':effect+=labels*delta
                    logs=np.log(z['base_mean'])+effect[:,None]*region[None]
                    mean=np.exp(logs-logsumexp(logs,axis=1,keepdims=True))*z['exposure'][:,None]
                    # Independent log-space normalization versus producer's ordinary sums.
                    budget=128*len(genes)*np.finfo(float).eps*np.maximum(1.,mean)
                    assert np.all(abs(mean-z['mean'])<=budget)
                    assert np.max(abs(z['mean'].sum(axis=1)-z['exposure']))<=128*len(genes)*np.finfo(float).eps*max(z['exposure'])
                    gamma=np.random.default_rng(np.random.SeedSequence([MASTER,phase,rep,3])).gamma(10.,.1,mean.shape)
                    counts=np.random.default_rng(np.random.SeedSequence([MASTER,phase,rep,4])).poisson(z['mean']*gamma)
                    assert np.array_equal(counts,X[query].toarray())
                    base=z['base_mean'].copy()
                signature=(X[:query[0]].copy(),base,ids.copy(),genes.copy())
                if rep in paired:
                    other=paired[rep]
                    assert (signature[0]!=other[0]).nnz==0 and all(np.array_equal(x,y) for x,y in zip(signature[1:],other[1:]))
                else:
                    regenerated=source_reference(rep,phase)
                    assert (signature[0]!=regenerated[0]).nnz==0
                    assert all(np.array_equal(x,y) for x,y in zip(signature[1:],regenerated[1:]))
                    paired[rep]=signature
                checked_cases.add(key)
            fit=case/arm/'fit';actual=independent_score(fit,truth[truth.specimen.eq('S1')])
            for key,value in actual.items():assert row[key]==value,(condition,rep,arm,key)
            assert file_record(row['fit_manifest']['path'])==row['fit_manifest']
            qc=json.loads((fit/'qc.json').read_text());assert qc['backend']=='cpu'
            trace=json.loads((case/arm/'trace.json').read_text())
            other='candidate' if arm=='comparator' else 'comparator'
            ot=json.loads((case/other/'trace.json').read_text())
            assert len(trace['gene'])==len(ot['gene'])==1
            with np.load(trace['gene'][0]['path']) as x,np.load(ot['gene'][0]['path']) as y:
                assert set(x.files)==set(y.files)
                for key in x.files:assert np.array_equal(x[key],y[key]),(condition,rep,key)
            audited.append(dict(condition=condition,rep=rep,arm=arm,**actual))
        summaries=[]
        for condition,arm in itertools.product(CONDITIONS,('comparator','candidate')):
            rows=[r for r in audited if r['condition']==condition and r['arm']==arm]
            upper=condition!='two_clone';event='null_error' if upper else 'exact_chr12_recovery'
            k=sum(r[event] for r in rows)
            result=dict(condition=condition,arm=arm,n=n,events=k,event=event,
                group_counts={str(g):sum(r['reported_groups']==g for r in rows) for g in sorted({r['reported_groups'] for r in rows})},
                mean_assignment_accuracy=sum(r['assignment_correct'] for r in rows)/(600*n),
                mean_tumour_call_fraction=sum(r['tumour_calls'] for r in rows)/(600*n))
            if not a.pilot:
                b,q=bound(k,n,upper,.05/6);result.update(simultaneous_bound=b,independent_binomial_bound=q,
                    direction='upper' if upper else 'lower',false_split_gate_passed=bool(b<=.05) if upper else None)
            summaries.append(result)
        write_json(dict(all_checks_passed=True,pilot=a.pilot,primary_evaluation_includes_pilot=False,
            specimen_conditions_checked=len(checked_cases),workflow_fits_checked=len(records),summaries=summaries,
            records=audited,biological_specimens=0,production_changed=False,scientific_criteria_changed=False,
            family=6,scope='Known-noise conditional simulation only; no biological/calibrated interval claim'),a.out/'report.json')
    except BaseException as error:
        finish_run(a.out,rec,error=error);raise
    finish_run(a.out,rec,backend='CPU independent full-workflow/generation/scoring audit')


if __name__=='__main__':main()
