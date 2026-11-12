"""Independent scalar generating means, target, operator and metric checks."""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose import simulate
from haplodose.provenance import start_run, finish_run, file_record


def average_ranks(values):
    order=sorted(range(len(values)),key=lambda i:values[i]);ranks=np.empty(len(values));start=0
    while start<len(order):
        end=start+1
        while end<len(order) and values[order[end]]==values[order[start]]:end+=1
        for index in order[start:end]:ranks[index]=(start+1+end)/2
        start=end
    return ranks


def correlation(x,y):
    n=len(x);a=math.fsum(x)/n;b=math.fsum(y)/n
    dx=[float(v)-a for v in x];dy=[float(v)-b for v in y]
    return math.fsum(v*w for v,w in zip(dx,dy))/math.sqrt(math.fsum(v*v for v in dx)*math.fsum(v*v for v in dy))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);panel=a.task/'panel_v1'
    inputs=dict(script=Path(__file__),manifest=panel/'run_manifest.json',
        values=panel/'primary_means_and_profiles.npz',regions=panel/'regions.json',
        metrics=panel/'metrics.tsv',report=panel/'report.json')
    previous=json.loads(inputs['manifest'].read_text());assert previous['status']=='complete'
    for group in ('input','output'):
        for record in previous[group]['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
    run=start_run(a.out,inputs,{},dict(scope='Scalar computational audit; no independent biological evidence'),mode='profile_pool_mechanism_audit')
    try:
        captured=[]
        def trace(frame,event,arg):
            if event=='return' and frame.f_code.co_filename==simulate.__file__ and frame.f_code.co_name=='draw':
                v=frame.f_locals
                # Recompute from rate, exposure and slope rather than using share.
                for L in v['L']:
                    unnormalized=[float(rate)*math.exp(math.log(float(L)/float(v['depth']))*float(g))
                                  for rate,g in zip(v['rate'],v['gamma'])]
                    total=math.fsum(unnormalized)
                    captured.append([float(L)*r/total for r in unnormalized])
        assert sys.getprofile() is None;sys.setprofile(trace)
        try:
            _,genes,meta,_,_,_=simulate.simulate(seed=20261109,n_chrom=4,genes_per_chrom=60,extra_genes=(),
                n_types=2,ref_cells_per_type=60,depth=1200.,specimens=[dict(name='A',n_cells=160,normal_frac=.1,
                clones=[dict(frac=.5,events=[('chr3',0.,1.,3)]),dict(frac=.5,events=[('chr3',0.,1.,4)])])])
        finally:sys.setprofile(None)
        M=np.array(captured)
        with np.load(inputs['values']) as z:v={k:z[k] for k in z.files}
        u=np.finfo(float).eps;operations=4*len(genes)+20;gamma=operations*u/(1-operations*u)
        error=np.abs(M-v['gene_mean'])
        assert (error<=gamma*np.maximum(np.abs(M),np.abs(v['gene_mean']))).all()
        source=a.task.parent/'profile_pool_export_20261109/parity_v1/fitted_new'
        used=pd.read_csv(source/'genes.tsv.gz',sep='\t');bins=pd.read_csv(source/'bins.tsv',sep='\t')
        columns=pd.Index(genes).get_indexer(used.gene);assert np.array_equal(columns,v['gene_columns'])
        mu=np.array([[math.fsum(float(row[columns[k]]) for k in np.flatnonzero(used.bin.eq(b)))
                      for b in bins.bin] for row in M])
        assert (np.abs(mu-v['used_gene_mean'])<=2*gamma*np.maximum(1,np.abs(mu))).all()
        query=v['query_mask'];ref=v['reference_mask']
        assert np.array_equal(ref,meta.role.eq('reference')&meta.ref_split.eq('fit')&meta.ref_group.eq('type0'))
        mean=np.array([math.fsum(mu[query,j]) for j in range(len(bins))])
        reference=np.array([math.fsum(mu[ref,j]) for j in range(len(bins))])
        truth=np.array([math.log(float(x)/math.fsum(mean)/(float(y)/math.fsum(reference))) for x,y in zip(mean,reference)])
        truth-=np.median(truth)
        assert (np.abs(truth-v['truth'])<=4*gamma*(1+np.abs(truth))).all()
        regions=json.loads(inputs['regions'].read_text());coverage=np.zeros(len(bins),int)
        for start,end in regions:
            assert start<end and bins.chromosome.iloc[start]==bins.chromosome.iloc[end-1]
            coverage[start:end]+=1
        assert np.array_equal(coverage,np.ones(len(bins),int))
        frame=pd.read_csv(inputs['metrics'],sep='\t',float_precision='round_trip');checked=0;max_difference=0.
        for row in frame.itertuples():
            profile=v[row.method+'__'+row.operator]
            if row.operator=='shared_naive_boundaries':
                original=v[row.method+'__none']
                for start,end in regions:
                    expected=math.fsum(original[start:end])/(end-start)
                    assert (np.abs(profile[start:end]-expected)<=gamma*(1+abs(expected))).all()
            # Stored target is independently checked above; use it here to isolate
            # score arithmetic from generator/target reference discrepancies.
            target=v['truth']
            scores=dict(MSE=math.fsum((float(x)-float(y))**2 for x,y in zip(profile,target))/len(target),
                Pearson=correlation(profile,target),Spearman=correlation(average_ranks(profile),average_ranks(target)))
            for name,value in scores.items():
                difference=abs(value-getattr(row,name));assert difference<=gamma*(1+abs(value))
                max_difference=max(max_difference,difference);checked+=1
        (a.out/'report.json').write_text(json.dumps(dict(generating_gene_means_checked=int(M.size),
            generating_bin_means_checked=int(mu.size),common_target_bins=len(truth),metric_checks=checked,
            scalar_max_metric_discrepancy=max_difference,scalar_max_mean_discrepancy=float(error.max()),
            regions=len(regions),scientific_criteria_unchanged=True,biological_validation=False),indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out,run,error=error);raise
    finish_run(a.out,run,backend='CPU scalar NumPy/math references')


if __name__=='__main__':main()
