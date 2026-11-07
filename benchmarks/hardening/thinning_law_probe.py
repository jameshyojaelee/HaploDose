"""Saved-proposal quadrature diagnostic for frozen-size/refitted-size thinning mismatch."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy import optimize,special,stats
from haplodose.clones import _fit_event_mixture
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--panel-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False);reports=[]
for seed in (20261108,20261109,20261110):
 root=a.root/f'trace_{seed}';events=json.loads((root/'trace.json').read_text());event=next(e for e in events if e['stage']=='event_mixtures');beta=next(e['beta'] for e in events if e['stage']=='difference_segments')
 with np.load(root/event['arrays']) as z:
  columns=z['cols'];mask=(columns>=209)&(columns<=217);O=z['Yp'][:,mask].copy();E=z['Ep'][:,mask].copy();R=np.broadcast_to(z['Rnb'],z['Yp'].shape)[:,mask].copy();r0=z['model_thin_size'][columns[mask]].copy();globalrows=z['global_cell_rows'][z['cells']].copy()
 nb=_fit_event_mixture(O,E,R,beta*np.log(1.5));truth=pd.read_csv(a.panel_root/f'seed_{seed}'/'data/cells.tsv',sep='\t').iloc[globalrows].truth_clone.to_numpy();tum=truth!='normal';target=truth=='clone2'
 def accuracy(labels):return max(float((labels[tum]==target[tum]).mean()),float((~labels[tum]==target[tum]).mean()))
 entry=dict(seed=seed,working_nb_gain=nb['gain'],working_nb_accuracy=accuracy(nb['responsibility']>.5),quadrature=[])
 for order in (32,64,128):
  nodes=[];weights=[]
  for size in r0:
   x,w=special.roots_jacobi(order,size/2-1,size/2-1);nodes.append((x+1)/2);weights.append(w/w.sum())
  nodes=np.array(nodes);weights=np.array(weights);r=2*R;fixed=special.gammaln(O+r)-special.gammaln(r)-special.gammaln(O+1)
  logM=np.log(2*E)[:,:,None]+np.log(nodes)[None];logr=np.log(r)[:,:,None]
  def component(eta):
   lm=logM+eta;den=np.logaddexp(logr,lm);terms=fixed[:,:,None]+r[:,:,None]*(logr-den)+O[:,:,None]*(lm-den)+np.log(weights)[None]
   perbin=special.logsumexp(terms,axis=2);q=np.exp(terms-perbin[:,:,None]);score=(q*(O[:,:,None]-(O+r)[:,:,None]*special.expit(lm-logr))).sum(axis=(1,2))
   return perbin.sum(axis=1),score
  null=optimize.minimize_scalar(lambda eta:-component(eta)[0].sum(),bracket=(nb['null_logtheta']-1,nb['null_logtheta']+1),options={'xtol':np.sqrt(np.finfo(float).eps)})
  baseline=component(null.x)[0]
  def objective(x):
   c0,s0=component(x[0]);c1,s1=component(x[1]);l0=c0-baseline-np.logaddexp(0,x[2]);l1=c1-baseline-np.logaddexp(0,-x[2]);resp=special.expit(l1-l0)
   return -np.logaddexp(l0,l1).sum()/len(O),np.array([-np.dot(1-resp,s0),-np.dot(resp,s1),np.sum(special.expit(x[2])-resp)])/len(O)
  fit=optimize.minimize(objective,nb['parameters'],jac=True,method='L-BFGS-B',bounds=[(None,None),(None,None),(-np.log(999),np.log(999))],options={'maxiter':1000,'ftol':np.finfo(float).eps,'gtol':np.finfo(float).eps**.75,'maxls':50})
  l0=component(fit.x[0])[0]-np.logaddexp(0,fit.x[2]);l1=component(fit.x[1])[0]-np.logaddexp(0,-fit.x[2]);labels=l1>l0
  entry['quadrature'].append(dict(order=order,optimizer_success=bool(fit.success),message=str(fit.message),gain=-float(fit.fun)*len(O),gradient=np.asarray(fit.jac).tolist(),accuracy=accuracy(labels),label_changes=int((labels!=(nb['responsibility']>.5)).sum()),parameters=fit.x.tolist()))
 reports.append(entry);print(json.dumps(entry),flush=True)
(a.out/'summary.json').write_text(json.dumps(dict(cases=reports,model='H|P ~ NB(mean=full_mean*P,size=1/refitted_phi), P~Beta(frozen_size/2,frozen_size/2); independent Gauss-Jacobi integration',scope='Saved proposal only. No model parameters or frozen outputs changed; truth used only for post-fit diagnostic accuracy. Nonconverged fits are reported unresolved, not accepted.'),indent=2))
