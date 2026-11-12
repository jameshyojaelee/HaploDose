"""Conditional binned discovery with physically withheld allele confirmation data."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import jax
import numpy as np
import pandas as pd

from haplodose.alleles import AlleleParams, SpecimenAlleles
from haplodose.clones import CloneParams, Rows, fit_specimen
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def array_hash(array):
    array=np.ascontiguousarray(array)
    return hashlib.sha256(str(array.dtype).encode()+json.dumps(array.shape).encode()+array.tobytes()).hexdigest()


def json_hash(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False).encode()).hexdigest()


def fit_and_freeze(y,a,d,snps,bins,selected,kappa,directory):
    directory.mkdir(exist_ok=False)
    params=CloneParams(n_programmes=0,spatial_neighbors=0,seed=20261108,split_mode='thin')
    allele=SpecimenAlleles(a[:,selected],d[:,selected],snps.iloc[selected],50.,AlleleParams(cell_kappa=kappa))
    allele.bins=bins
    reset=np.zeros(40,bool);reset[0]=True
    result=fit_specimen(Rows(y,np.arange(120)),y.sum(axis=1,dtype=np.float64),
                        np.full((40,1),1/40),reset,np.full(40,.15),params,alleles=allele)
    model=result['model']
    arrays={name:np.array(result[name],copy=True) for name in ('R','LL','doses','pi','Yk','Lk','XLk')}
    for name in ('post','path','dose','E'):
        arrays['profile_'+name]=np.array(result['prof'][name],copy=True)
    for name in ('Ak','Dk','B','allele_post','phi','offset','gamma','U'):
        value=getattr(model,name)
        if value is not None:
            arrays[name]=np.array(value,copy=True)
    for value in arrays.values():
        value.setflags(write=False)
    regions={}
    for k in range(len(result['doses'])):
        for run in model.alleles.runs(model.allele_post[k],model.Dk[k]):
            if run['resolved']:
                key=(run['chromosome'],run['start'],run['end'])
                regions[key]=dict(chromosome=key[0],start=key[1],end=key[2])
    meta=dict(selected_snp_indices=selected.tolist(),reference_kappa=50.,
              fitted_kappa=float(allele.kappa),
              fitted_nu=float(allele.nu) if allele.nu is not None and np.isfinite(allele.nu) else None,
              regions=list(regions.values()),evidence=result['evidence'])
    write_json(model.log,directory/'search_log.json')
    meta['search_log_sha256']=file_record(directory/'search_log.json')['sha256']
    write_json(meta,directory/'discovery.json')
    np.savez_compressed(directory/'discovery.npz',**arrays)
    fingerprint=dict(arrays={name:array_hash(value) for name,value in arrays.items()},
                     metadata=json_hash(meta))
    write_json(fingerprint,directory/'fingerprint.json')
    return arrays,meta,fingerprint


def reserved_summaries(arrays,meta,a,d,snps,reserved,directory):
    """Read reserved counts only after candidate/region/labels are frozen.

    Raw summaries only: no adaptive threshold, nuisance fit, test significance,
    new group, or overwrite of discovery arrays.
    """
    labels=arrays['R'].argmax(axis=1)
    eligible=[k for k in range(arrays['R'].shape[1]) if np.sum(labels==k)>=30]
    results=[]
    for region in meta['regions']:
        positions=snps.pos.to_numpy();chrom=snps.chromosome.to_numpy()
        sites=reserved[(chrom[reserved]==region['chromosome'])&(positions[reserved]>region['start'])
                       &(positions[reserved]<=region['end'])]
        for left,right in itertools.combinations(eligible,2):
            aa=np.array([a[labels==k][:,sites].sum(axis=0) for k in (left,right)])
            dd=np.array([d[labels==k][:,sites].sum(axis=0) for k in (left,right)])
            results.append(dict(region=region,candidate_pair=[left,right],snps=sites,
                                n_cells=[int(np.sum(labels==k)) for k in (left,right)],
                                alternate_counts=aa,total_counts=dd,
                                status='raw_reserved_summary_only'))
    write_json(dict(contrasts=results,reason='' if results else 'no_frozen_region_and_eligible_pair',
                    inference='No p-value, calibration or biological-validation claim'),directory/'reserved.json')


def descriptive_refit(arrays,a,d,snps,bins,kappa,directory):
    # Fresh object prevents full-data nu/kappa updates from mutating discovery.
    full=SpecimenAlleles(a,d,snps,50.,AlleleParams(cell_kappa=kappa))
    full.bins=bins
    ak,dk=full.pool(arrays['R'])
    post,b=full.fit(ak,dk)
    np.savez_compressed(directory/'descriptive_all_snp.npz',post=post,B=b,Ak=ak,Dk=dk)
    write_json(dict(scope='Descriptive fixed-membership all-SNP refit',kappa=float(full.kappa),
                    nu=float(full.nu) if np.isfinite(full.nu) else None),directory/'descriptive.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--cell-model',choices=('binomial','beta'),required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    assert jax.default_backend()=='cpu'
    kappa=None if args.cell_model=='binomial' else 2.
    run=start_run(args.out,dict(script=Path(__file__),design=args.task/'DESIGN.md'),
                  dict(cell_model=args.cell_model,cells=120,bins=40,snps=40,seed=20261108),
                  dict(scope='Complete fit_specimen discovery at fixed binned inputs/reference; not public raw-input or biological validation'),
                  mode='allele_discovery_workflow')
    try:
        snps=pd.DataFrame(dict(chromosome='chr1',pos=np.arange(40)*1000+1,cM=np.arange(40)*.001))
        bins=pd.DataFrame(dict(chromosome='chr1',start=np.arange(40)*1000,end=(np.arange(40)+1)*1000))
        rng=np.random.default_rng([20261108,810])
        effect=np.zeros(40);effect[:10]=.6*np.log(1.5)
        rate=np.exp(effect);rate/=rate.sum()
        y=rng.negative_binomial(1/.15,1/(1+.15*800*rate),size=(120,40)).astype(np.float32)
        d=np.full((120,40),20.)
        a=np.repeat(np.repeat([10.,12.,16.],40)[:,None],40,axis=1)
        proposal=np.arange(0,40,2);reserved=np.arange(1,40,2)
        np.savez_compressed(args.out/'fixed_inputs.npz',Y=y,A_even=a[:,proposal],D=d,
                            proposal_indices=proposal,reserved_indices=reserved)
        scenarios=[('odd_balanced',10.,False),('odd_same',19.,False),('odd_reverse',1.,False),
                   ('absent_even_balanced',10.,True),('absent_even_strong',19.,True)]
        records=[];first={}
        for name,odd,absent in scenarios:
            ac=a.copy();dc=d.copy();ac[:,reserved]=odd
            if absent:
                ac[:,proposal]=0.;dc[:,proposal]=0.
            np.savez_compressed(args.out/(name+'_alleles.npz'),A=ac,D=dc)
            arms=('discovery_subset',) if absent else ('current_all_snps','discovery_subset')
            for arm in arms:
                print(json.dumps(dict(scenario=name,arm=arm,status='starting')),flush=True)
                directory=args.out/(name+'_'+arm)
                selected=proposal if arm=='discovery_subset' else np.arange(40)
                arrays,meta,fingerprint=fit_and_freeze(y,ac,dc,snps,bins,selected,kappa,directory)
                protected={p.name:file_record(p)['sha256'] for p in directory.iterdir()}
                reserved_summaries(arrays,meta,ac,dc,snps,reserved,directory)
                descriptive_refit(arrays,ac,dc,snps,bins,kappa,directory)
                assert all(file_record(directory/name)['sha256']==digest for name,digest in protected.items())
                assert {key:array_hash(value) for key,value in arrays.items()}==fingerprint['arrays']
                group_sizes=np.bincount(arrays['R'].argmax(axis=1),minlength=arrays['R'].shape[1])
                record=dict(scenario=name,arm=arm,discovery=fingerprint,
                            counts_sha256=array_hash(y),selected_A_sha256=array_hash(ac[:,selected]),
                            selected_D_sha256=array_hash(dc[:,selected]),group_sizes=group_sizes,
                            frozen_after_confirmation_and_refit=True)
                records.append(record)
                key=(arm,absent)
                if key not in first:
                    first[key]=fingerprint
                elif arm=='discovery_subset':
                    assert fingerprint==first[key], 'Reserved SNPs changed complete discovery'
        write_json(records,args.out/'cases.json')
        # Repeat the baseline after all unrelated scenarios; never carry fitted
        # full-data concentration, memberships or phase state between calls.
        repeat_a=a.copy();repeat_a[:,reserved]=10.
        _,_,repeat=fit_and_freeze(y,repeat_a,d,snps,bins,proposal,kappa,args.out/'baseline_repeat')
        assert repeat==first[('discovery_subset',False)], 'Independent fit order changed discovery'
        write_json(dict(discovery_fits=len(records)+1,cell_model=args.cell_model,
                        complete_binned_search_invariance=True,fit_order_invariant=True,
                        descriptive_refit_isolated=True,biological_specimens=0,
                        production_changed=False,calibrated=False),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
