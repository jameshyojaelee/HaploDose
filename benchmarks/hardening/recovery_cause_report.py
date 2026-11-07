"""Saved-decision audit and independent latent-information counterexample; compute only."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import integrate,special,stats
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
reports=[]
for seed in (20261108,20261109,20261110):
 r=a.root/f'trace_{seed}'; events=json.loads((r/'trace.json').read_text());bins=pd.read_csv(r/'fit/bins.tsv',sep='\t');region=bins.loc[bins.chromosome.eq('chr12'),'bin'].to_numpy();lo=int(region[0]);hi=int(region[len(region)//2-1]);key=f'_bins{lo}-{hi}'
 c=pd.read_csv(r/'assertion_view/cells.tsv.gz',sep='\t');normal=c[c.specimen.eq('S1') & c.truth_clone.eq('normal')]
 normal.to_csv(a.out/f'normal_cells_{seed}.tsv',sep='\t',index=False)
 selected=[]
 for e in events:
  ctx=str(e['context'])
  if e['stage']=='event_mixtures':
   selected.append(dict(stage=e['stage'],index=e['index'],found=[x for x in e['found'] if key in x['source']],returned=e['returned']))
  elif key in ctx:
   item={k:v for k,v in e.items() if k!='arrays'}
   if e['stage']=='difference_segments':
    with np.load(r/e['arrays']) as z:
     keep=(z['idx']>=lo)&(z['idx']<=hi)
     item['beta']=e['beta']; item['step']=float(e['beta']*np.log(1.5)); item['maximum_nonzero_posterior']=float(z['post'][:,keep][:,:,[0,2]].max()); item['adjusted_mean']=float(z['d'][:,keep].mean()); item['mean_variance']=float(z['var'][:,keep].mean()); item['target_bins']=z['idx'][keep].tolist();item['target_d']=z['d'][:,keep].tolist();item['target_variance']=z['var'][:,keep].tolist();item['target_posterior']=z['post'][:,keep].tolist()
   selected.append(item)
 reports.append(dict(seed=seed,gate=json.loads((r/'summary.json').read_text())['gate'],normal_false_calls=int(normal.malignant_call.sum()),normal_n=len(normal),chr12_first_half=[lo,hi],events=selected))
(a.out/'causes.json').write_text(json.dumps(reports,indent=2))
# Original effect/noise, no model fitting or cutoff selection. A balanced equal-variance
# Gaussian mixture has density derivative proportional to -x+a*tanh(a*x/sigma²).
# a<sigma implies that derivative is strictly negative for x>0: one continuous mode.
delta=.6*np.log(1.5);sigma=.25;half=delta/2
v=sigma*sigma+half*half
logp=lambda x: special.logsumexp([stats.norm.logpdf(x,-half,sigma),stats.norm.logpdf(x,half,sigma)])-np.log(2)
logq=lambda x: stats.norm.logpdf(x,0,np.sqrt(v))
results=[]
for tolerance in (1e-8,1e-10,1e-12):
 kl,err=integrate.quad(lambda x: np.exp(logp(x))*(logp(x)-logq(x)),-np.inf,np.inf,epsabs=tolerance,epsrel=tolerance)
 results.append(dict(requested_precision=tolerance,KL=kl,quadrature_error=err,pinsker_power_upper_at_alpha05=min(1.,.05+np.sqrt(600*(kl+err)/2))))
# Reinterpret the same regional latent effect as transcription, keeping normalized
# rates and NB emissions equal. This is an algebraic counterexample, not new truth.
rng=np.random.default_rng(20261108);state=rng.normal(0,sigma,600);labels=np.r_[np.zeros(300),np.ones(300)]
genomic=state+delta*labels;transcription=genomic.copy()
(a.out/'identification.json').write_text(json.dumps(dict(delta=delta,state_sd=sigma,unimodality_derivative_upper=-1+(half/sigma)**2,quadrature=results,identical_regional_effect=np.array_equal(genomic,transcription),scope='KL bound applies to IID latent mixture versus moment-matched Gaussian, not the original fixed-balanced generator. Exact state/CNA reinterpretation preserves original fixed-balanced emissions. No biological replication, threshold change or new evaluation claim.'),indent=2))
