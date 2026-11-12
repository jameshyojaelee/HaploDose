"""Independent scalar fixtures and all-case membership-confidence recount."""
import argparse
from fractions import Fraction
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.provenance import start_run, finish_run, file_record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);panel=a.task/'panel_v1'
    source=Path(__file__).with_name('clone_membership_confidence.py')
    spec=importlib.util.spec_from_file_location('confidence_fixture',source)
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    inputs=dict(script=Path(__file__),helper=source,manifest=panel/'run_manifest.json')
    previous=json.loads(inputs['manifest'].read_text());assert previous['status']=='complete'
    for record in previous['output']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
    for name in ('cases.tsv','reliability.tsv','assignments.json','summary.tsv','report.json'):inputs[name]=panel/name
    run=start_run(a.out,inputs,{},dict(scope='Independent computational recount, no biological validation'),
                  mode='clone_membership_confidence_audit')
    try:
        R=np.array([[4,4,0],[7,0,1],[0,6,2],[1,2,5]],float)/8
        i,j,q,_=helper.coassignment(R)
        exact=np.array([float(sum(Fraction(int(R[x,k]*8),8)*Fraction(int(R[y,k]*8),8)
                                for k in range(3))) for x,y in zip(i,j)])
        assert np.array_equal(q,exact)
        permutation_checks=0
        reference=R@R.T
        for cells,components in itertools.product(itertools.permutations(range(4)),itertools.permutations(range(3))):
            cells=np.array(cells);components=np.array(components)
            x,y,actual,_=helper.coassignment(R[cells][:,components])
            assert np.array_equal(actual,reference[cells[x],cells[y]])
            permutation_checks+=1
        rejected=0
        for invalid in (np.array([[np.nan,0],[1,0]]),np.array([[np.inf,0],[1,0]]),
                        np.array([[-.1,1.1],[1,0]]),np.array([[0.,0],[1,0]]),np.array([[1.,1],[1,0]])):
            try:helper.coassignment(invalid)
            except ValueError:rejected+=1
            else:raise AssertionError('Invalid responsibilities accepted')
        assert rejected==5
        for q0,same,expected in ((np.array([1.,0.]),np.array([True,False]),(0.,1.,0.)),
                                 (np.array([0.,1.]),np.array([True,False]),(1.,1.,1.)),
                                 (np.array([.5,.5]),np.array([True,False]),(.25,0.,0.))):
            m=helper.metrics(q0,same)
            assert (m['brier'],m['confident_coverage'],m['confident_unconditional_error'])==expected
        frame=pd.read_csv(inputs['cases.tsv'],sep='\t',float_precision='round_trip')
        reliability=pd.read_csv(inputs['reliability.tsv'],sep='\t',float_precision='round_trip')
        inventory=json.loads(inputs['assignments.json'].read_text())
        assert len(frame)==len(inventory)==1200 and len(reliability)==12000
        assert not frame[['arm','role','replicate']].duplicated().any()
        i,j=np.triu_indices(600,1);checks=0;maximum_brier_difference=0.;scalar_checks=0
        for index,record in enumerate(inventory):
            assert file_record(record['path'])['sha256']==record['sha256']
            key=(record['arm'],record['role'],record['replicate'])
            row=frame.loc[frame.arm.eq(key[0])&frame.role.eq(key[1])&frame.replicate.eq(key[2])].iloc[0]
            with np.load(record['path']) as z:raw=z['raw_R'];truth=z['simulated_CNA_labels']
            R=raw.astype(float);R/=R.sum(axis=1)[:,None]
            # Matrix multiplication instead of the panel's pairwise einsum.
            q=(R@R.T)[i,j];same=truth[i]==truth[j]
            n=len(q);K=R.shape[1];u=np.finfo(float).eps
            gamma=(8*K+16)*u/(1-(8*K+16)*u)
            brier=float(np.sum(np.square(q-same))/n)
            difference=abs(brier-row.brier);assert difference<=gamma
            maximum_brier_difference=max(maximum_brier_difference,difference)
            if record['replicate'] in (0,199):
                scalar=math.fsum((float(v)-int(t))**2 for v,t in zip(q,same))/n
                assert abs(scalar-row.brier)<=gamma;scalar_checks+=1
            high=q>=.95;low=q<=.05;selected=high|low
            wrong=(high&~same)|(low&same)
            counts=dict(pairs=n,confident_pairs=int(selected.sum()),confident_errors=int(wrong.sum()),
                high_confidence_separations=int(low.sum()),high_confidence_same=int(high.sum()))
            for name,value in counts.items():assert value==row[name],(key,name)
            for name,value in dict(confident_coverage=selected.sum()/n,
                confident_unconditional_error=wrong.sum()/n,
                confident_unconditional_correct=(selected.sum()-wrong.sum())/n,
                one_population_brier=float(np.mean(~same))).items():assert abs(value-row[name])<=gamma
            rel=reliability.loc[reliability.arm.eq(key[0])&reliability.role.eq(key[1])&reliability.replicate.eq(key[2])]
            decile=np.minimum((q*10).astype(int),9)
            for d in range(10):
                r=rel.loc[rel.decile.eq(d)].iloc[0];mask=decile==d
                assert int(mask.sum())==r.pairs and int(same[mask].sum())==r.same_pairs
                assert abs(float(q[mask].sum())-r.probability_sum)<=gamma*max(1,int(mask.sum()))
            checks+=1
        summary=pd.read_csv(inputs['summary.tsv'],sep='\t',float_precision='round_trip')
        for _,row in summary.iterrows():
            f=frame.loc[frame.arm.eq(row.arm)&frame.role.eq(row.role)]
            assert len(f)==row.specimens==200
            assert int(f.confident_pairs.sum())==row.total_confident_pairs
            assert int(f.confident_errors.sum())==row.total_confident_errors
            assert int(f.high_confidence_separations.gt(0).sum())==row.specimens_with_confident_separation
        (a.out/'report.json').write_text(json.dumps(dict(cases_recounted=checks,
            pair_reliability_strata_recounted=12000,exact_fraction_pairs=len(exact),
            permutation_checks=permutation_checks,invalid_inputs_rejected=rejected,
            scalar_full_pair_Brier_checks=scalar_checks,max_Brier_discrepancy=maximum_brier_difference,
            scientific_criteria_unchanged=True,biological_validation=False),indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out,run,error=error);raise
    finish_run(a.out,run,backend='CPU independent NumPy/Fraction/scalar references')


if __name__=='__main__':main()
