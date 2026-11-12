"""Independently audit merge eligibility, latent-profile identity and frozen parity."""
import argparse
import itertools
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from haplodose import clones
from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    source=a.task/'trace_v1';selected=json.loads((source/'report.json').read_text())
    inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md',selected=source/'report.json')
    prior=a.task.parent/'gene_veto_fresh_20261109'
    for path in sorted(source.rglob('*')):
        if path.is_file() and path.suffix in ('.npz','.json','.tsv','.gz','.npy'):
            inputs['trace:'+str(path.relative_to(source))]=path
    original=[]
    for folder in sorted(prior.glob('chunk_*_v1')):
        inputs['original:'+folder.name]=folder/'results.json';original.extend(json.loads((folder/'results.json').read_text()))
    for row in original:
        case=Path(row['case']);key=f"{row['condition']}:{row['rep']}:{row['arm']}"
        for name in ('clones.tsv','cells.tsv.gz','differences.tsv'):
            inputs[key+':'+name]=case/row['arm']/'fit'/name
    for row in selected:
        for name in ('posteriors.npy','profiles.tsv.gz','segments.tsv'):
            inputs['prior:'+row['condition']+':'+str(row['rep'])+':'+name]=Path(row['prior'])/'candidate/fit'/name
    rec=start_run(a.out,inputs,{},dict(scope='Selected mechanisms and descriptive frozen-panel decomposition; no new event definition'),mode='gene_state_trace_audit')
    try:
        cases=[]
        for row in selected:
            out=Path(row['out']);old=Path(row['prior'])/'candidate/fit'
            assert row['all_original_arrays_and_tables_identical'] and file_record(row['trace']['path'])==row['trace']
            for name in ('cells.tsv.gz','clones.tsv','differences.tsv','segments.tsv','profiles.tsv.gz'):
                pd.testing.assert_frame_equal(pd.read_csv(out/'fit'/name,sep='\t'),pd.read_csv(old/name,sep='\t'),check_exact=True)
            assert np.array_equal(np.load(out/'fit/posteriors.npy'),np.load(old/'posteriors.npy'),equal_nan=True)
            trace=json.loads((out/'trace.json').read_text());matched=[x for x in trace['weighted'] if x['stage']=='merge_membership'];mi=0;merges=[]
            for entry in trace['stages']:
                z=np.load(entry['arrays'])
                if entry['stage']=='agglomerate':
                    assert len(entry['members'])==len(z['Y_after'])
                    assert sorted(sum(entry['members'],[]))==list(range(len(z['Y_before'])))
                    for i,m in enumerate(entry['members']):assert np.array_equal(z['Y_before'][m].sum(axis=0),z['Y_after'][i])
                elif entry['stage']=='merge_assignment':
                    R=z['proposal_R'];labels=R.argmax(axis=1);counts=np.bincount(labels,minlength=R.shape[1])
                    eligible=[k for k in range(1,R.shape[1]) if counts[k]>=entry['min_cells']]
                    pairs=[list(x) for x in itertools.combinations(eligible,2)]
                    if pairs:
                        m=matched[mi];mi+=1
                        assert m['pairs']==pairs and np.array_equal(m['labels'],labels)
                        assert np.array_equal(m['proposal_responsibility'],R) and np.array_equal(m['full_responsibility'],z['full_R'])
                    merges.append(dict(index=entry['index'],components=R.shape[1],proposal_MAP_counts=counts.tolist(),
                        full_MAP_counts=np.bincount(z['full_R'].argmax(axis=1),minlength=R.shape[1]).tolist(),
                        full_soft_mass=z['full_R'].sum(axis=0).tolist(),eligible=eligible,compared_pairs=pairs))
            assert mi==len(matched)
            final=next(x for x in trace['stages'] if x['stage']=='final');z=np.load(final['arrays']);R=z['R'];doses=z['doses']
            counts=np.bincount(R.argmax(axis=1),minlength=R.shape[1]);order=[0]+sorted(range(1,R.shape[1]),key=lambda k:-counts[k])
            P=np.load(out/'fit/posteriors.npy');assert np.array_equal(P[z['rows'],:R.shape[1]],R[:,order].astype(np.float16))
            table=pd.read_csv(out/'fit/clones.tsv',sep='\t').sort_values('column')
            assert table.n_cells_map.tolist()==counts[order].tolist()
            budget=8*len(R)*np.finfo(R.dtype).eps*np.maximum(1.,R.sum(axis=0)[order])
            assert np.all(abs(table.n_eff.to_numpy()-R.sum(axis=0)[order])<=budget)
            distances=[dict(a=i,b=j,exactly_identical=bool(np.array_equal(doses[i],doses[j])),
                max_absolute_difference=float(np.max(abs(doses[i]-doses[j])))) for i,j in itertools.combinations(range(1,len(doses)),2)]
            cases.append(dict(condition=row['condition'],rep=row['rep'],all_frozen_exports_identical=True,
                reported_tumour_profiles=len(doses)-1,final_MAP_counts=counts.tolist(),final_soft_mass=R.sum(axis=0).tolist(),
                dose_pairs=distances,merge_calls=merges,initial_agglomerations=[dict(context=e['context'],members=e['members'])
                    for e in trace['stages'] if e['stage']=='agglomerate'],
                terminal_evidence=final['evidence'],model_log=final['log']))
        decomposition=[]
        for row in original:
            fit=Path(row['case'])/row['arm']/'fit'
            table=pd.read_csv(fit/'clones.tsv',sep='\t');tum=table[~table.normal_like]
            cells=pd.read_csv(fit/'cells.tsv.gz',sep='\t');cells=cells[cells.specimen.eq('S1')]
            assert int(table.n_cells_map.sum())==len(cells)==600
            for item in table.itertuples():assert int(item.n_cells_map)==sum(cells.clone==item.clone)
            assert len(tum)==row['reported_groups']
            diff=pd.read_csv(fit/'differences.tsv',sep='\t');diff=diff[diff.clone_a.isin(tum.clone)&diff.clone_b.isin(tum.clone)]
            decomposition.append(dict(condition=row['condition'],rep=row['rep'],arm=row['arm'],reported=len(tum),
                original_event=row['null_error'],nonempty_MAP_profiles=int(sum(tum.n_cells_map>0)),
                MAP_at_least30_profiles=int(sum(tum.n_cells_map>=30)),zero_MAP_profiles=int(sum(tum.n_cells_map==0)),
                exported_tumour_differences=len(diff)))
        # Minimal original behavior: four identical tumour components survive
        # because only one has proposal MAP support. No scientific tolerance.
        m=clones.SpecimenModel.__new__(clones.SpecimenModel);m.thin=True;m.params=SimpleNamespace(min_cells=30);m.log=[];m.Ak=None
        likelihood=np.tile(np.array([-100.,0.,0.,0.,0.]),(600,1));m.loglik=lambda doses,fold:likelihood
        m.prop_fold='prop'
        pi=np.array([.01,.396,.297,.198,.099]);R=np.tile(np.array([0.,.4,.3,.2,.1]),(600,1))
        dose=np.r_[np.full(5,.2),np.zeros(5)];doses=[np.zeros(10)]+[dose.copy() for _ in range(4)]
        returned=m.merge_unsupported(R,doses,pi)
        assert not returned[2] and len(returned[0])==5
        assert all(np.array_equal(returned[0][1],d) for d in returned[0][2:])
        counterexample=dict(tumour_components_before=4,tumour_components_after=4,merged=False,
            all_tumour_profiles_exactly_identical=True,interpretation='Original merge can leave exactly nonidentifiable components; not proof every zero-MAP component is invalid')
        write_json(dict(all_checks_passed=True,cases=cases,decomposition=decomposition,exact_counterexample=counterexample,
            original_state_candidate_errors=41,original_state_evaluation_n=100,
            scientific_endpoint_changed=False,production_changed=False,biological_specimens=0),a.out/'report.json')
    except BaseException as error:finish_run(a.out,rec,error=error);raise
    finish_run(a.out,rec,backend='CPU independent identities, pair eligibility, exact profile counterexample and frozen census')


if __name__=='__main__':main()
