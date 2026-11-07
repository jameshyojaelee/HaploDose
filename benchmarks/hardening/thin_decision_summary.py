"""Summarize saved chr12 HMM decisions without rerunning inference; compute only."""
import argparse
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from haplodose.io import new_output_dir, write_json
from haplodose.provenance import start_run, finish_run


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--trace-root',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    root=args.trace_root.resolve()
    trace=json.loads((root/'trace.json').read_text())
    bins=pd.read_csv(root/'fit/bins.tsv',sep='\t').set_index('bin')
    region=bins.index[bins.chromosome.eq('chr12')].to_numpy()
    lo,hi=int(region[0]),int(region[len(region)//2-1])
    selected=[]
    for entry in trace:
        if entry['stage']!='difference_segments':
            continue
        source=entry.get('context',{}).get('source','')
        match=re.search(r'_bins(\d+)-(\d+)$',source)
        if match and tuple(map(int,match.groups()))==(lo,hi):
            selected.append(entry)
    if not selected:
        raise ValueError('No saved first-half chr12 difference decision; no replacement selected')
    out=new_output_dir(args.out)
    inputs=dict(trace=root/'trace.json',trace_manifest=root/'run_manifest.json',
                bins=root/'fit/bins.tsv',runner=Path(__file__))
    for entry in selected:
        inputs[f'capture_{entry["index"]}']=root/entry['arrays']
    event=next((e for e in trace if e['stage']=='event_mixtures'),None)
    if event:
        inputs['saved_mixture_state']=root/event['arrays']
    record=start_run(out,inputs,{},dict(purpose='Saved path/posterior decision audit; no HMM rerun'),mode='thin_decision_summary')
    try:
        reports=[]
        for entry in selected:
            with np.load(root/entry['arrays']) as z:
                idx=z['idx'];reset=z['sub_reset'].astype(bool)
                for pair in range(len(z['path'])):
                    signed=np.array([-1,0,1])[z['path'][pair]]
                    starts=np.flatnonzero(np.r_[True,(signed[1:]!=signed[:-1])|reset[1:]])
                    for a,b in zip(starts,np.r_[starts[1:],len(signed)]):
                        full=idx[a:b]
                        if not np.any((full>=lo)&(full<=hi)):
                            continue
                        direction=int(signed[a]);support=float(z['post'][pair,a:b,direction+1].mean())
                        d=z['d'][pair,a:b];var=z['var'][pair,a:b];weight=1/var
                        mean=float((weight*d).sum()/weight.sum())
                        params=entry['parameters']
                        state_posterior=z['post'][pair,a:b,direction+1]
                        credible=state_posterior>=params['posterior']
                        core_starts=np.flatnonzero(credible & np.r_[True,~credible[:-1]])
                        core_ends=np.flatnonzero(credible & np.r_[~credible[1:],True])+1
                        cores=[]
                        for ca,cb in zip(core_starts,core_ends):
                            core_bins=full[ca:cb]
                            core_weight=1/var[ca:cb]
                            core_mean=float((core_weight*d[ca:cb]).sum()/core_weight.sum())
                            cores.append(dict(first_bin=int(core_bins[0]),last_bin=int(core_bins[-1]),
                                start=int(bins.loc[int(core_bins[0]),'start']),
                                end=int(bins.loc[int(core_bins[-1]),'end']),n_bins=len(core_bins),
                                passes_existing_min_bins=len(core_bins)>=params['min_bins'],
                                nonneutral=direction!=0,
                                qualifies=direction!=0 and len(core_bins)>=params['min_bins'],
                                minimum_state_posterior=float(state_posterior[ca:cb].min()),
                                mean_state_posterior=float(state_posterior[ca:cb].mean()),
                                mean_adjusted_logratio=core_mean,
                                z=float(core_mean*np.sqrt(core_weight.sum()))))
                        reports.append(dict(capture=entry['index'],context=entry['context'],pair=pair,
                            full_first_bin=int(full[0]),full_last_bin=int(full[-1]),
                            chromosome=str(bins.loc[int(full[0]),'chromosome']),
                            start=int(bins.loc[int(full[0]),'start']),end=int(bins.loc[int(full[-1]),'end']),
                            signed_state=direction,run_bins=len(full),support=support,
                            per_bin_state_posterior=[dict(bin=int(bin_id),posterior=float(prob))
                                for bin_id,prob in zip(full,state_posterior)],
                            maximal_contiguous_credible_cores=cores,
                            required_bins=params['min_bins'],required_posterior=params['posterior'],
                            length_pass=len(full)>=params['min_bins'],posterior_pass=support>=params['posterior'],
                            nonneutral_pass=direction!=0,mean_adjusted_logratio=mean,
                            mean_variance=float(var.mean()),tau=float(z['tau'][pair]),
                            beta=entry['beta'],step=float(entry['beta']*np.log(1.5)),
                            z=float(mean*np.sqrt(weight.sum())),original_results=entry['results'][pair],
                            parameters=params))
        variance=None
        if event:
            with np.load(root/event['arrays']) as z:
                phi=z['model_phi'];r0=z['model_thin_size'];s=.5
                actual=phi+((1-s)/s)*(1+phi)/(r0+1)
                variance=dict(trace_event=event['index'],model_phi_range=[float(phi.min()),float(phi.max())],
                    frozen_size_range=[float(r0.min()),float(r0.max())],
                    max_abs_model_phi_minus_inverse_size=float(abs(phi-1/r0).max()),
                    max_abs_implied_half_phi_minus_assumed=float(abs(actual-phi/s).max()),
                    derivation='For NB Y with mean m,var=m+phi*m² and independent P~Beta(s*r0,(1-s)*r0), H|Y,P~Binomial(Y,P): Var(H)=s*m+s²*phi*m²+s*(1-s)*(1+phi)*m²/(r0+1). Hence phi_half=phi+((1-s)/s)*(1+phi)/(r0+1), equal phi/s iff r0=1/phi.',
                    limit='Model phi is fitted, not known true dispersion; no correction applied.')
        write_json(dict(region_bins=[lo,hi],runs=reports,variance=variance,
            interpretation='Reported HMM states/posteriors are saved outputs. No inference rerun and no criteria changed.'),out/'summary.json')
        finish_run(out,record,backend='cpu')
    except Exception as error:
        finish_run(out,record,backend='cpu',error=error)
        raise


if __name__=='__main__':
    main()
