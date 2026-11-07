"""Independent bounded NB-mixture diagnosis from saved thinning arrays, on compute."""
import argparse
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy.optimize import minimize, minimize_scalar
from scipy.special import expit, gammaln

from haplodose.io import new_output_dir, write_json
from haplodose.provenance import start_run, finish_run


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--trace-root', required=True, type=Path)
    ap.add_argument('--candidate-source', help='Exact saved event source; default first half of chr12')
    ap.add_argument('--out', required=True, type=Path)
    args = ap.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    root = args.trace_root.resolve()
    trace = json.loads((root / 'trace.json').read_text())
    bins = pd.read_csv(root / 'fit/bins.tsv', sep='\t')
    chr12 = bins.loc[bins.chromosome.eq('chr12'), 'bin'].to_numpy()
    first, last = int(chr12[0]), int(chr12[len(chr12) // 2 - 1])
    selected = None
    for entry in trace:
        if entry['stage'] != 'event_mixtures':
            continue
        for index, candidate in enumerate(entry['found']):
            match = re.search(r'_bins(\d+)-(\d+)$', candidate['source'])
            wanted = candidate['source'] == args.candidate_source if args.candidate_source else (
                match is not None and (int(match[1]), int(match[2])) == (first, last))
            if wanted:
                selected = entry, index, candidate
                break
        if selected:
            break
    if selected is None:
        raise ValueError('Requested saved candidate not present; no replacement region selected')
    entry, index, candidate = selected
    path = root / entry['arrays']
    out = new_output_dir(args.out)
    rec = start_run(out, dict(runner=Path(__file__), trace=root / 'trace.json',
        trace_manifest=root / 'run_manifest.json', arrays=path, bins=root / 'fit/bins.tsv'),
        dict(candidate_source=candidate['source'], trace_event=entry['index']),
        dict(purpose='Independent optimization of saved proposal only; no re-inference or truth selection'),
        mode='thin_saved_mixture_optimization')
    try:
        with np.load(path) as z:
            lo, hi = map(int, re.search(r'_bins(\d+)-(\d+)$', candidate['source']).groups())
            mask = (z['cols'] >= lo) & (z['cols'] <= hi)
            O, E = z['Yp'][:, mask].copy(), z['Ep'][:, mask].copy()
            R = np.broadcast_to(z['Rnb'], z['Yp'].shape)[:, mask].copy()
            original_labels = z[f'found_groups_{index}'].copy()
            phi, r0 = z['model_phi'].copy(), z['model_thin_size'].copy()
        if not np.isfinite(E).all() or np.any(E <= 0) or np.any(R <= 0):
            raise ValueError('Saved expectations/sizes must be positive finite')
        constant = (gammaln(O + R) - gammaln(R) - gammaln(O + 1)).sum(axis=1)
        n = len(O)

        def component(logtheta):
            mu = E * np.exp(logtheta)
            ll = (O * np.log(mu / (mu + R)) + R * np.log(R / (mu + R))).sum(axis=1)
            score = (R * (O - mu) / (R + mu)).sum(axis=1)
            return ll, score

        def objective(x, offset=0.):
            a, da = component(x[0]); b, db = component(x[1])
            pi = expit(x[2])
            la = a - np.logaddexp(0., x[2]); lb = b - np.logaddexp(0., -x[2])
            total = np.logaddexp(la, lb)
            rb = np.exp(lb - total)
            gradient = np.array([((1-rb)*da).sum(), (rb*db).sum(), (rb-pi).sum()])
            return -(total.sum() + offset) / n, -gradient / n

        def rate(labels):
            return minimize_scalar(lambda x: -component(x)[0][labels].sum(), bounds=(-8.,8.),
                                   method='bounded', options={'xatol': 1e-12}).x
        null = rate(np.ones(n, bool)); null_ll = component(null)[0].sum()
        starts = [np.array([rate(original_labels == 0), rate(original_labels == 1),
                           np.log(max((original_labels == 1).mean(), 1e-6) / max((original_labels == 0).mean(),1e-6))])]
        starts += [np.array([null-d, null+d, p]) for d,p in [(0.05,0.),(0.2,0.),(0.5,0.),(0.2,-1.),(0.2,1.)]]
        fits = []
        for offset_name, offset in [('omitted_constants',0.), ('full_NB_constants',float(constant.sum()))]:
            for start in starts:
                fit = minimize(lambda x: objective(x, offset), start, jac=True, method='L-BFGS-B',
                    bounds=[(-8.,8.),(-8.,8.),(-12.,12.)],
                    options={'maxiter':3000,'ftol':1e-14,'gtol':1e-9,'maxls':50})
                x = fit.x
                a,_=component(x[0]);b,_=component(x[1]);pi=expit(x[2])
                rb=np.exp(b-np.logaddexp(0.,-x[2])-np.logaddexp(a-np.logaddexp(0.,x[2]),b-np.logaddexp(0.,-x[2])))
                gain = -objective(x)[0]*n-null_ll
                fits.append(dict(offset=offset_name,start=start.tolist(),parameters=x.tolist(),
                    theta=np.exp(x[:2]).tolist(),pi=float(pi),gain=float(gain),success=bool(fit.success),
                    message=str(fit.message),iterations=int(fit.nit),score_per_cell=objective(x)[1].tolist(),
                    max_abs_score_per_cell=float(np.abs(objective(x)[1]).max()),
                    agreement_original=max(float(np.mean((rb>.5)==original_labels)),float(np.mean((rb<=.5)==original_labels))),
                    at_bounds=bool(np.any(np.isclose(np.abs(x),[8.,8.,12.],atol=1e-5)))))
        best=max(fits,key=lambda x:x['gain'])
        x=np.asarray(best['parameters'])
        finite_difference=[]
        for axis in range(3):
            delta=np.zeros(3);delta[axis]=1e-5
            finite_difference.append(float((objective(x+delta)[0]-objective(x-delta)[0])/(2e-5)))
        a,_=component(x[0]);b,_=component(x[1])
        la=a-np.logaddexp(0.,x[2]);lb=b-np.logaddexp(0.,-x[2])
        posterior=np.exp(lb-np.logaddexp(la,lb))
        np.savez_compressed(out/'independent_region.npz',O=O,E=E,Rnb=R,
            original_labels=original_labels,posterior=posterior,labels=(posterior>.5),parameters=x)
        s=.5
        phi_half=phi+((1-s)/s)*(1+phi)/(r0+1)
        np.savez_compressed(out/'variance.npz',model_phi=phi,frozen_thin_size=r0,
                            implied_half_phi=phi_half,assumed_half_phi=phi/s)
        write_json(dict(candidate=candidate,trace_event=entry['index'],shape=list(O.shape),fits=fits,best=best,
            original_gain=float(candidate['gain']),gain_improvement=float(best['gain']-candidate['gain']),
            derivative_check=dict(analytic=objective(x)[1].tolist(),central_difference=finite_difference,
                                  central_difference_step=1e-5),
            original_stop_thresholds=dict(omitted=1e-6*max(abs(null_ll+candidate['gain']),1.),
                with_full_NB_constants=1e-6*max(abs(null_ll+candidate['gain']+constant.sum()),1.)),
            derivation=['Y has mean m and variance m+phi*m^2; P independent Beta(s*r0,(1-s)*r0).',
                'E[P]=s, Var(P)=s*(1-s)/(r0+1); H|Y,P is Binomial(Y,P).',
                'Var(H)=s*m+s^2*phi*m^2+s*(1-s)*(1+phi)*m^2/(r0+1).',
                'phi_half=phi+((1-s)/s)*(1+phi)/(r0+1); equals phi/s only when r0=1/phi.'],
            variance_summary=dict(max_abs_phi_difference=float(np.max(abs(phi_half-phi/s))),
                max_abs_phi_minus_inverse_frozen_size=float(np.max(abs(phi-1/r0)))),
            interpretation='Optimization and variance diagnostics only. Variance formula assumes the underlying Y has the stated NB phi; fitted model_phi is an estimate. Multi-start solutions are not a proof of global optimum. No acceptance criteria changed, no correction applied, no truth used.'),out/'summary.json')
        finish_run(out,rec,backend='cpu')
    except Exception as error:
        finish_run(out,rec,backend='cpu',error=error)
        raise


if __name__=='__main__':
    main()
