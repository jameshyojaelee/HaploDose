"""Independent conditional expectation of hard and margin-selected RNA pools."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from haplodose.provenance import start_run,finish_run

p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--trace-root',type=Path,required=True)
a=p.parse_args();a.out.mkdir(exist_ok=False)
rec=start_run(a.out,{'script':__file__},{},{'scope':'Conditional fitted-probability expectation, not measured purity or calibrated uncertainty'},mode='weight_expectation')
rows=[]
try:
 for seed in (20261108,20261109,20261110):
    root=a.trace_root/f'trace_{seed}'
    for event in json.loads((root/'trace.json').read_text()):
        if event['stage']!='event_fit':continue
        z=np.load(root/event['arrays']);q=z['q'];g=q>.5;exposure=z['E'].sum(axis=1);theta=np.exp(z['parameters'][:2])
        values=[]
        for mode,weights in [('hard',np.ones(len(q))),('margin',np.abs(2*q-1))]:
            purity=[];mean=[];se=[];support=[]
            for group in (False,True):
                take=g==group;w=weights[take];E=exposure[take];prob=q[take]
                purity.append(float((w*E*prob).sum()/(w*E).sum()))
                rate=((1-prob)*theta[0]+prob*theta[1])
                mean.append(float((w*E*rate).sum()/(w*E).sum()))
                support.append(float(w.sum()**2/(w*w).sum()))
            direct=float(np.log(mean[1]/mean[0]))
            alternate=float(logsumexp(z['parameters'][:2],b=[1-purity[1],purity[1]])-
                            logsumexp(z['parameters'][:2],b=[1-purity[0],purity[0]]))
            # Two independent formulas on the saved fit, tolerance bounded by
            # a few dozen floating operations at IEEE double precision.
            np.testing.assert_allclose(direct,alternate,rtol=64*np.finfo(float).eps,atol=64*np.finfo(float).eps)
            values.append(dict(mode=mode,purity=purity,predicted_log_contrast=direct,effective_cells=support,
                               independent_formula_discrepancy=alternate-direct))
        rows.append(dict(seed=seed,trace_index=event['index'],latent_contrast=float(np.diff(z['parameters'][:2])[0]),pools=values))
 (a.out/'summary.json').write_text(json.dumps(rows,indent=2));finish_run(a.out,rec,backend='cpu')
except Exception as error:
 finish_run(a.out,rec,error=error);raise
