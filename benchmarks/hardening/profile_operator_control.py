"""Exact shared-operator RNA profile control; boundaries never use DNA truth."""
import argparse,hashlib,itertools,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm,spearmanr
from haplodose.profile import ProfileParams

CASES={'2313287_STOMACH':('pseudobulk','hd_gene_ratio'),
       '42MGBA_CENTRAL_NERVOUS_SYSTEM':('pseudobulk','hd_gene_ratio'),
       'CRC.SW.U0001.T':('naive','hd_ratio'),'CRC.SW.U0008.T':('naive','hd_ratio'),
       'TCGA-05-4389-01':('naive','hd_lr')}

def partition(y,penalty,min_size):
    y=np.asarray(y,float);n=len(y)
    if n<min_size:return np.full(n,y.mean()),[(0,n)],float(((y-y.mean())**2).sum()+penalty)
    centered=y-y.mean();s=np.r_[0.,np.cumsum(centered)];ss=np.r_[0.,np.cumsum(centered**2)]
    cost=np.full(n+1,np.inf);previous=np.full(n+1,-1,int);cost[0]=0
    for stop in range(min_size,n+1):
        start=np.arange(stop-min_size+1);size=stop-start
        loss=np.maximum(0.,ss[stop]-ss[start]-(s[stop]-s[start])**2/size)
        choices=cost[start]+loss+penalty;j=int(np.argmin(choices))
        cost[stop]=choices[j];previous[stop]=start[j]
    intervals=[];end=n;fitted=np.empty(n)
    while end:
        begin=previous[end]
        if begin<0:raise ValueError('No valid partition')
        intervals.append((int(begin),int(end)));fitted[begin:end]=y[begin:end].mean();end=begin
    return fitted,intervals[::-1],float(cost[n])

def verify():
    y=np.array([.1,-.2,.3,1.,1.1,.8,.9]);penalty=.17
    for minimum in (1,2,3):
        ref=[]
        for bits in itertools.product((False,True),repeat=len(y)-1):
            ends=[0]+[i+1 for i,b in enumerate(bits) if b]+[len(y)]
            if min(np.diff(ends))<minimum:continue
            ref.append(sum(float(((y[a:b]-y[a:b].mean())**2).sum())+penalty for a,b in zip(ends[:-1],ends[1:])))
        actual=partition(y,penalty,minimum)[2]
        np.testing.assert_allclose(actual,min(ref),rtol=64*np.finfo(float).eps,atol=64*np.finfo(float).eps)
        fitted,regions,_=partition(y,penalty,minimum)
        scale=np.log(2.)
        converted,converted_regions,_=partition(scale*y,scale**2*penalty,minimum)
        assert regions==converted_regions
        np.testing.assert_allclose(converted,scale*fitted,rtol=64*np.finfo(float).eps,atol=64*np.finfo(float).eps)

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);verify();results=[]
    for case,(naive,hd) in CASES.items():
        source=a.inputs/f'{case}_profiles.tsv.gz';frame=pd.read_csv(source,sep='\t')
        if frame.bin.duplicated().any() or not frame.bin.is_monotonic_increasing:raise ValueError('Unordered/ambiguous bins')
        # subarm.type_windows uses log2(CPM+1); Kinker and U-CAN use ln.
        # Correlation is invariant to this conversion, shared penalties are not.
        units={naive:'log2' if case=='TCGA-05-4389-01' else 'ln',hd:'ln'}
        values={method:frame[method].to_numpy()*(np.log(2.) if units[method]=='log2' else 1.) for method in (naive,hd)}
        good=np.isfinite(frame[naive]) & np.isfinite(frame[hd]);idx=np.flatnonzero(good)
        breaks=np.r_[True,(np.diff(idx)!=1)|(frame.arm.to_numpy()[idx[1:]]!=frame.arm.to_numpy()[idx[:-1]])|
                     (np.diff(frame.bin.to_numpy()[idx])!=1)]
        blocks=np.split(idx,np.flatnonzero(breaks)[1:])
        if not len(idx):raise ValueError('No common RNA support')
        diffs=np.concatenate([np.diff(values[naive][b]) for b in blocks if len(b)>1])
        sigma=np.median(np.abs(diffs-np.median(diffs)))/(np.sqrt(2)*norm.ppf(.75))
        if not np.isfinite(sigma) or sigma<=0:raise ValueError('Noise scale unidentifiable; no epsilon substitution')
        penalty=3*sigma**2*np.log(len(idx));minimum=ProfileParams().min_bins
        boundaries={}
        for method in (naive,hd):
            output=np.full(len(frame),np.nan);boundaries[method]=[]
            for b in blocks:
                fit,regions,_=partition(values[method][b],penalty,minimum);output[b]=fit
                boundaries[method]+=[(int(frame.bin.iloc[b[x]]),int(frame.bin.iloc[b[y-1]])) for x,y in regions]
            frame[method+'_own_partition']=output
        # Truth enters only after every operator decision has been made.
        common=good & frame.common_callable & np.isfinite(frame.truth)
        scores={}
        for method in (naive,hd,naive+'_own_partition',hd+'_own_partition'):
            v=frame.loc[common,method];truth=frame.loc[common,'truth']
            scores[method]=dict(n=len(v),spearman=float(spearmanr(v,truth).statistic))
        record=dict(case=case,input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),methods=[naive,hd],
                    input_log_units=units,partition_output_units='ln',
                    sigma_from_naive_RNA=float(sigma),penalty=float(penalty),minimum_segment_bins=minimum,
                    boundaries=boundaries,scores=scores,biological_n=1,
                    scope='Locked development diagnostic; Gaussian RNA log-ratio working loss, not count uncertainty calibration or independent biological validation')
        frame.to_csv(a.out/f'{case}_profiles.tsv.gz',sep='\t',index=False);results.append(record)
    (a.out/'summary.json').write_text(json.dumps(results,indent=2));print(json.dumps(results),flush=True)

if __name__=='__main__':main()
