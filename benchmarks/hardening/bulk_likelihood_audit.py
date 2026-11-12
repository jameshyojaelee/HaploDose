"""Independent scalar count sums, SciPy PLN integrals and HMM recursion audit."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import logsumexp, softmax

from haplodose.provenance import start_run, finish_run, file_record


def log_integral(y, lm, sd, relative):
    """Mode-scaled infinite-domain integral, omitting shared log(y!) constant.

    The positive-count Poisson part omits y log(y)-y as well. Both constants
    cancel in state likelihood ratios. Root finding and adaptive integration
    are independent of the production Newton/Hermite implementation.
    """
    inv = 1/(sd*sd)
    def slope(o): return y-np.exp(min(lm+o,700.))-o*inv
    lo,hi = -max(1.,sd),max(1.,sd)
    while slope(lo)<0: lo *= 2
    while slope(hi)>0: hi *= 2
    mode = brentq(slope,lo,hi,xtol=np.finfo(float).tiny,rtol=4*np.finfo(float).eps)
    rate = np.exp(lm+mode); width = 1/np.sqrt(rate+inv)
    if y:
        u = lm+mode-np.log(y)
        top = -y*(np.expm1(u)-u)-mode*mode*inv/2
    else: top = -rate-mode*mode*inv/2
    def integrand(z):
        delta = z*width
        if delta>700: return 0.
        # Taylor remainder avoids cancellation at the tightly concentrated mode.
        remainder = np.expm1(delta)-delta
        log_change = (y-rate-mode*inv)*delta-rate*remainder-delta*delta*inv/2
        return np.exp(log_change) if log_change>-745 else 0.
    value,error = quad(integrand,-np.inf,np.inf,epsabs=0.,epsrel=relative,limit=400)
    return top+np.log(width*value)-np.log(sd),error/value


def independent_hmm(emit, reset):
    """Separate float64 probability transitions and chromosome-wise recursion."""
    prior = np.array([.2,.5,.2,.1])
    transition = np.empty((4,4))
    for a in range(4):
        for b in range(4): transition[a,b] = 1-.002 if a==b else .002*prior[b]/(1-prior[a])
    t = np.log(transition); out = np.empty_like(emit,dtype=float)
    for block in np.split(np.arange(len(emit)),np.flatnonzero(reset)[1:]):
        e = emit[block].astype(float); n=len(block)
        a=np.zeros_like(e); b=np.zeros_like(e)
        a[0]=np.log(prior)+e[0]; a[0]-=logsumexp(a[0])
        for k in range(1,n):
            a[k]=e[k]+logsumexp(a[k-1,:,None]+t,axis=0); a[k]-=logsumexp(a[k])
        for k in range(n-2,-1,-1):
            b[k]=logsumexp(t+e[k+1,None,:]+b[k+1,None,:],axis=1); b[k]-=max(b[k])
        out[block]=softmax(a+b,axis=1)
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    args=p.parse_args();args.out.mkdir(exist_ok=False)
    root=args.task.parent.parent; trace=args.task/'trace_v1'
    hd=root/'dev/bulk_posterior_export_20261109/stage_v1/haplodose'
    raw=root/'dev/hardening_20261106_2119/profile_v6/raw_refit_000/counts.tsv'
    inputs=dict(script=Path(__file__),design=args.task/'AUDIT_DESIGN.md',manifest=trace/'run_manifest.json',
        windows=trace/'loss_windows.tsv',mapping=trace/'window_gene_bins.tsv',counts=raw,
        bins=hd/'bins.tsv',genes=hd/'genes.tsv.gz',samples=hd/'samples.tsv',posterior=hd/'rna_state_posterior.npy')
    run=start_run(args.out,inputs,dict(quadrature_relative_targets=[1e-8,1e-11]),
        dict(scope='Independent numerical audit of exposed-panel mechanism trace'),mode='bulk_likelihood_audit')
    try:
        manifest=json.loads((trace/'run_manifest.json').read_text());assert manifest['status']=='complete'
        for record in manifest['output']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        windows=pd.read_csv(trace/'loss_windows.tsv',sep='\t',float_precision='round_trip')
        mapping=pd.read_csv(trace/'window_gene_bins.tsv',sep='\t')
        bins=pd.read_csv(hd/'bins.tsv',sep='\t');genes=pd.read_csv(hd/'genes.tsv.gz',sep='\t')
        meta=pd.read_csv(hd/'samples.tsv',sep='\t');counts=pd.read_csv(raw,sep='\t',index_col=0)
        saved=np.load(hd/'rna_state_posterior.npy');reset=np.r_[True,bins.chromosome.to_numpy()[1:]!=bins.chromosome.to_numpy()[:-1]]
        # Declared centromere/gap blocks are independent chains too.
        reset[1:] |= bins.genomic_block.to_numpy()[1:]!=bins.genomic_block.to_numpy()[:-1]
        numerics=[];hmm=[];replay={};count_checks=0
        for sample,groups in mapping.groupby('sample',sort=True):
            with np.load(trace/(sample+'.npz')) as z: fit={k:z[k] for k in z.files}
            expected=[]
            for _,group in genes.groupby('bin',sort=True):expected.append(sum(int(counts.at[g,sample]) for g in group.gene))
            assert np.array_equal(expected,fit['Y']);count_checks+=len(expected)
            ref=independent_hmm(fit['emit'],reset);i=meta['sample'].tolist().index(sample)
            discrepancy=float(np.max(np.abs(ref-saved[i])))
            # Existing HMM emission-offset test criterion is retained, not relaxed.
            assert discrepancy<1e-4, (sample,discrepancy)
            hmm.append(dict(sample=sample,bins=len(ref),max_absolute_posterior_discrepancy=discrepancy))
            local=np.empty((len(ref),4));local[:]=np.nan
            for b in sorted(groups.bin.unique()):
                y=float(np.float32(fit['Y'][b]));sd=float(np.float32(fit['sd'][b]))
                coarse=[];fine=[];errors=[]
                for state in range(4):
                    lm=float(fit['means'][b,state])
                    a,ea=log_integral(y,lm,sd,1e-8);c,ec=log_integral(y,lm,sd,1e-11)
                    coarse.append(a);fine.append(c);errors.append(max(ea,ec)+abs(a-c))
                fine=np.array(fine);coarse=np.array(coarse)
                production=fit['emit'][b].astype(float)
                # Subtract neutral; factorial and other shared constants cancel.
                delta=(production-production[1])-(fine-fine[1])
                local[b]=softmax(fine+np.log([.2,.5,.2,.1]))
                numerics.append(dict(sample=sample,bin=int(b),count=y,sd=sd,
                    max_state_log_bf_discrepancy=float(np.max(np.abs(delta))),
                    max_quadrature_log_error_estimate=max(errors),
                    max_refinement_log_bf_change=float(np.max(np.abs((fine-fine[1])-(coarse-coarse[1])))),
                    loss_neutral_log_bf_reference=float(fine[0]-fine[1]),
                    loss_neutral_log_bf_production=float(production[0]-production[1])))
            replay[sample]=(ref,local,fit)
        pd.DataFrame(numerics).to_csv(args.out/'quadrature_bins.tsv',sep='\t',index=False)
        pd.DataFrame(hmm).to_csv(args.out/'hmm_replay.tsv',sep='\t',index=False)
        records=[]
        for row in windows.itertuples():
            group=mapping[mapping['sample'].eq(row.sample)&mapping.window.eq(row.window)]
            bi=group.bin.to_numpy(int);ref,local,fit=replay[row.sample]
            production=fit['local'][bi,0]
            error=float(np.max(np.abs(production-local[bi,0])))
            ref_median=float(np.median(local[bi,0]));smooth_median=float(np.median(ref[bi,0]))
            records.append(dict(sample=row.sample,window=row.window,category=row.category,
                local_loss_reference=ref_median,smoothed_loss_reference=smooth_median,
                local_probability_max_discrepancy=error,
                local_above_half=bool(ref_median>.5),smoothed_above_half=bool(smooth_median>.5),
                beta_loss_amplitude=float(fit['beta']*-np.log(2)),posterior_call_possible_from_state_range=bool(
                    fit['beta']*-np.log(2)<-.356426864862442)))
        checked=pd.DataFrame(records);checked.to_csv(args.out/'window_reference.tsv',sep='\t',index=False)
        groups={}
        for category,frame in checked.groupby('category'):
            groups[category]=dict(windows=len(frame),local_above_half=int(frame.local_above_half.sum()),
                smoothed_above_half=int(frame.smoothed_above_half.sum()),
                local_above_smoothed_below=int((frame.local_above_half&~frame.smoothed_above_half).sum()),
                score_threshold_below_attainable_state_range=int((~frame.posterior_call_possible_from_state_range).sum()))
        n=pd.DataFrame(numerics)
        report=dict(scalar_count_checks=count_checks,patients=len(replay),loss_windows=len(checked),
            quadrature_bins=len(n),quadrature_integrals=8*len(n),
            maximum_log_bf_discrepancy=float(n.max_state_log_bf_discrepancy.max()),
            maximum_refinement_change=float(n.max_refinement_log_bf_change.max()),
            maximum_quad_log_error_estimate=float(n.max_quadrature_log_error_estimate.max()),
            maximum_local_probability_discrepancy=float(checked.local_probability_max_discrepancy.max()),
            maximum_hmm_probability_discrepancy=max(r['max_absolute_posterior_discrepancy'] for r in hmm),
            categories=groups,probability_half_is_diagnostic_not_a_call_threshold=True,
            quadrature_error_estimates_are_not_formal_bounds=True,method_promotion=False,biological_calibration=False)
        (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        for record in run['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
