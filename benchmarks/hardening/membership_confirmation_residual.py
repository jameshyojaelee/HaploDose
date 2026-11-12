"""Does the saved proposal mixture predict its confirmation counts? Diagnostic only."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import nbinom
from haplodose.provenance import start_run,finish_run

p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
p.add_argument('--trace-root',type=Path,required=True);p.add_argument('--panel-root',type=Path,required=True)
a=p.parse_args();a.out.mkdir(exist_ok=False)
rec=start_run(a.out,{'script':__file__},{},{'scope':'Selected-parent conditional diagnostic; no independent uncertainty or biological calibration'},mode='confirmation_residual')
report=[]
try:
 for seed in (20261108,20261109,20261110):
    root=a.trace_root/f'trace_{seed}'
    truth=pd.read_csv(a.panel_root/f'seed_{seed}'/'data/cells.tsv',sep='\t')
    for event in json.loads((root/'trace.json').read_text()):
        if event['stage']!='event_fit':continue
        z=np.load(root/event['arrays']);q=z['q'];E=z['E'];r=np.broadcast_to(z['r'],E.shape)
        theta=np.exp(z['parameters'][:2]);mu=np.stack([E*x for x in theta]);prob=np.stack([1-q,q])
        means=mu.sum(axis=2)
        within=(mu+mu**2/r).sum(axis=2)
        variance=(prob*within).sum(axis=0)+q*(1-q)*(means[1]-means[0])**2
        # Independent SciPy NB moments and total-variance decomposition.
        sm,sv=nbinom.stats(r[None],r[None]/(r[None]+mu),moments='mv')
        marginal=(prob*sm.sum(axis=2)).sum(axis=0)
        reference=(prob*(sv.sum(axis=2)+(sm.sum(axis=2)-marginal[None])**2)).sum(axis=0)
        np.testing.assert_allclose(variance,reference,rtol=64*np.finfo(float).eps,atol=64*np.finfo(float).eps)
        actual=truth.iloc[z['rows'][z['cells']]].truth_clone.to_numpy()
        tumour=np.isin(actual,['clone1','clone2']);y=(actual=='clone2').astype(float)
        row=dict(seed=seed,trace_index=event['index'],selected_cells=len(q),
                 selected_normal_cells=int((actual=='normal').sum()),
                 brier_tumour=float(np.mean((q[tumour]-y[tumour])**2)),
                 interpretation='Postfit truth diagnostic only; q is fitted on selected parents, not calibrated DNA membership',pools=[])
        for mode,w in [('hard',np.ones(len(q))),('margin',np.abs(2*q-1))]:
            groups=[]
            for group in (0,1):
                use=(q>.5)==bool(group);ww=w[use];exposure=E[use].sum(axis=1)
                expected=float(ww@marginal[use]);observed=float(ww@z['Ytest'][use].sum(axis=1))
                v=float(ww**2@variance[use]);denominator=float(ww@exposure)
                groups.append(dict(group=group,observed_weighted_total=observed,expected_weighted_total=expected,
                    conditional_sd=float(np.sqrt(v)),conditional_standardized_residual=(observed-expected)/np.sqrt(v),
                    observed_rate=observed/denominator,predicted_rate=expected/denominator,
                    fitted_purity=float((ww*exposure*q[use]).sum()/denominator),
                    postfit_truth_purity=float((ww*exposure*y[use]).sum()/denominator)))
            row['pools'].append(dict(mode=mode,groups=groups,
                observed_log_contrast=float(np.log(groups[1]['observed_rate']/groups[0]['observed_rate'])),
                predicted_log_contrast=float(np.log(groups[1]['predicted_rate']/groups[0]['predicted_rate']))))
        report.append(row)
 (a.out/'summary.json').write_text(json.dumps(report,indent=2));finish_run(a.out,rec,backend='cpu')
except Exception as error:
 finish_run(a.out,rec,error=error);raise
