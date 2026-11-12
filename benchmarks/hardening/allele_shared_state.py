"""Falsify genomic interpretation of SNP holdout under shared cell-level ASE."""
import argparse
import importlib.util
import itertools
import json
from pathlib import Path

import jax
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import linear_sum_assignment

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def sibling(name):
    spec=importlib.util.spec_from_file_location(name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def reserved_test(arrays,meta,a,d,snps,reserved,rng):
    labels=arrays['R'].argmax(axis=1)
    eligible=[k for k in range(arrays['R'].shape[1]) if np.sum(labels==k)>=30]
    pairs=list(itertools.combinations(eligible,2));regions=[]
    for region in meta['regions']:
        mask=(snps.chromosome.to_numpy()[reserved]==region['chromosome'])
        pos=snps.pos.to_numpy()[reserved]
        sites=reserved[mask&(pos>region['start'])&(pos<=region['end'])]
        if len(sites):
            regions.append(dict(region=region,sites=sites))
    # All depths are equal and positive in this fixed diagnostic.
    fractions=np.column_stack([a[:,r['sites']].sum(axis=1)/d[:,r['sites']].sum(axis=1)
                               for r in regions]) if regions else np.empty((120,0))
    permutations=np.vstack([np.arange(120)]+[rng.permutation(120) for _ in range(999)])
    maxima=np.zeros(1000)
    if pairs and regions:
        for left,right in pairs:
            means=[fractions[permutations[:,labels==k]].mean(axis=1) for k in (left,right)]
            maxima=np.maximum(maxima,np.max(np.abs(means[0]-means[1]),axis=1))
    available=bool(pairs and regions)
    p=float((1+np.sum(maxima[1:]>=maxima[0]))/1000) if available else 1.
    return dict(available=available,p=p,observed=float(maxima[0]),eligible=eligible,
                regions=regions,pairs=pairs),dict(labels=labels,fractions=fractions,
                permutations=permutations,maxima=maxima)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--cell-model',choices=('binomial','beta'),required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    assert jax.default_backend()=='cpu'
    workflow=sibling('allele_discovery_workflow');roles=sibling('allele_site_roles')
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',
                workflow=Path(workflow.__file__),roles=Path(roles.__file__))
    run=start_run(args.out,inputs,dict(replicates=40,cell_model=args.cell_model,permutations=999),
                  dict(scope='Conditional SNP-holdout mechanism falsification; not biological validation'),
                  mode='allele_shared_state')
    try:
        snps=pd.DataFrame(dict(chromosome='chr1',pos=np.arange(40)*1000+1,cM=np.arange(40)*.001))
        bins=pd.DataFrame(dict(chromosome='chr1',start=np.arange(40)*1000,end=(np.arange(40)+1)*1000))
        selected=np.array([s for s,row in snps.iterrows() if roles.site_role(row.chromosome,row.pos,genome_build='hg38',seed=20261108)=='discovery'])
        reserved=np.setdiff1d(np.arange(40),selected)
        write_json(dict(discovery=selected,reserved=reserved),args.out/'roles.json')
        kappa=None if args.cell_model=='binomial' else 2.
        records=[]
        for rep in range(40):
            streams=[np.random.default_rng(s) for s in np.random.SeedSequence([20261108,911,rep]).spawn(4)]
            rna_rng,state_rng,count_rng,permutation_rng=streams
            effect=np.zeros(40);effect[:10]=.6*np.log(1.5);rate=np.exp(effect);rate/=rate.sum()
            y=rna_rng.negative_binomial(1/.15,1/(1+.15*800*rate),size=(120,40)).astype(np.float32)
            truth=state_rng.permutation(np.repeat([0,1],60))
            site_state=state_rng.integers(0,2,size=(120,40))
            # Common uniforms couple conditions without changing their marginal laws.
            uniforms=count_rng.random((120,40,20))
            probs=dict(constant=np.full((120,40),.65),site_independent=.5+.3*site_state,
                       shared_state=np.broadcast_to((.5+.3*truth)[:,None],(120,40)))
            for condition,prob in probs.items():
                directory=args.out/f'{rep:03d}_{condition}'
                a=(uniforms<prob[:,:,None]).sum(axis=2).astype(float);d=np.full((120,40),20.)
                arrays,meta,fingerprint=workflow.fit_and_freeze(y,a,d,snps,bins,selected,kappa,directory)
                protected={p.name:file_record(p)['sha256'] for p in directory.iterdir()}
                test,reference=reserved_test(arrays,meta,a,d,snps,reserved,permutation_rng)
                np.savez_compressed(directory/'audit_inputs.npz',Y=y,A=a,D=d,truth=truth,**reference)
                write_json(test,directory/'reserved_test.json')
                assert all(file_record(directory/name)['sha256']==sha for name,sha in protected.items())
                labels=reference['labels'];table=np.array([[np.sum((labels==k)&(truth==t)) for t in (0,1)] for k in range(arrays['R'].shape[1])])
                rows,cols=linear_sum_assignment(-table)
                records.append(dict(replicate=rep,condition=condition,available=test['available'],p=test['p'],
                                    positive=test['p']<=.05,group_sizes=np.bincount(labels,minlength=arrays['R'].shape[1]),
                                    latent_all_cell_matching=float(table[rows,cols].sum()/120),
                                    frozen_discovery=fingerprint))
            print(json.dumps(dict(replicate=rep,status='complete')),flush=True)
        summaries=[]
        for condition in probs:
            subset=[r for r in records if r['condition']==condition];hits=sum(r['positive'] for r in subset);n=len(subset)
            summaries.append(dict(condition=condition,n=n,positives=hits,available=sum(r['available'] for r in subset),
                simultaneous_lower=float(stats.beta.ppf(.05/6,hits,n-hits+1)) if hits else 0.,
                simultaneous_upper=float(stats.beta.ppf(1-.05/6,hits+1,n-hits)) if hits<n else 1.,
                mean_latent_matching=float(np.mean([r['latent_all_cell_matching'] for r in subset]))))
        write_json(records,args.out/'cases.json')
        write_json(dict(cell_model=args.cell_model,summaries=summaries,biological_specimens=0,
                        scope='Candidate held-out contrasts; not production whole-specimen split decisions',
                        genomic_vs_ASE_identifiable=False,production_changed=False),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
