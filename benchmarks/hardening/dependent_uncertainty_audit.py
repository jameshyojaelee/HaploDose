"""Independent calibration-rank, interval identity and coverage/power recount."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize,stats


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True);p.add_argument('--previous',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);directory=a.task/'panel_v2'
    manifest=json.loads((directory/'run_manifest.json').read_text())
    if manifest['status']!='complete':raise AssertionError('Incomplete experiment')
    for r in manifest['input']['files'].values():
        if hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()!=r['sha256']:raise AssertionError('Input changed')
    summary=json.loads((directory/'summary.json').read_text())
    if any(r['log_p']>0 or r['log_p_holm']>0 for r in summary['coverage']):
        raise AssertionError('Invalid log probability')
    # The second execution repairs only probability reporting. Scientific
    # intervals, bank draws and the six endpoint decisions must be unchanged.
    first_run=a.task/'panel_v1'
    for name in ('intervals.tsv','calibration.json'):
        if (directory/name).read_bytes()!=(first_run/name).read_bytes():
            raise AssertionError('Reporting-only correction changed scientific artifacts')
    original=pd.read_csv(a.previous/'intervals.tsv',sep='\t',float_precision='round_trip')
    frame=pd.read_csv(directory/'intervals.tsv',sep='\t',float_precision='round_trip')
    keys=['phase','condition','effect','interval','replicate']
    old=original[original.method.eq('simultaneous')].sort_values(keys).reset_index(drop=True)
    retained=frame[frame.method.eq('simultaneous')].sort_values(keys).reset_index(drop=True)
    pd.testing.assert_frame_equal(old,retained,check_dtype=False,check_exact=True)
    new=frame[frame.method.eq('dependent_maximum')].sort_values(keys).reset_index(drop=True)
    if len(new)!=10800:raise AssertionError('Missing candidate intervals')
    for column in keys+['first_bin','last_bin','true_mu_I','estimated_tau','estimated_beta']:
        pd.testing.assert_series_equal(new[column],old[column],check_dtype=False,check_exact=True)
    if not np.isfinite(new[['lower','upper']]).all().all() or (new.lower>new.upper).any():
        raise AssertionError('Invalid candidate intervals in this positive-count panel')
    np.testing.assert_array_equal(new.covered,(new.lower<=new.true_mu_I)&(new.true_mu_I<=new.upper))
    np.testing.assert_array_equal(new.excludes_zero,(new.lower>0)|(new.upper<0))
    banks=np.load(directory/'calibration_maxima.npz');simulations=np.load(a.previous/'simulations.npz')
    calibration=json.loads((directory/'calibration.json').read_text());ranks=[]
    for condition,r in calibration.items():
        values=banks[condition];n=len(values);k=r['rank'];delta=.0125/3;q=1-.0125
        if n!=20000 or k!=1+int(stats.binom.ppf(1-delta,n,q)):raise AssertionError('Independent tolerance rank differs')
        failure=float(stats.beta.cdf(q,k,n+1-k));previous=float(stats.beta.cdf(q,k-1,n+2-k))
        if failure>delta or previous<=delta:raise AssertionError('Nonminimal/invalid tolerance rank')
        if np.sort(values)[k-1]!=r['critical']:raise AssertionError('Calibration cutoff differs from bank')
        # Reconstruct three preselected Gaussian draws. Evaluate every eligible
        # interval with fsum, without prefix arithmetic or calibration helper.
        S=simulations[f'coverage_{condition}_effect0_covariance'];L=np.linalg.cholesky(S)
        rng=np.random.default_rng(r['seed']);draws=rng.standard_normal((20000,120))
        discrepancies=[]
        for index in (0,10000,19999):
            Z=draws[index]@L.T;reference=0.
            for block in range(0,120,30):
                for first in range(block,block+30):
                    for stop in range(first+1,block+31):
                        variance=math.fsum(S[first:stop,first:stop].ravel())
                        reference=max(reference,abs(math.fsum(Z[first:stop]))/math.sqrt(variance))
            difference=float(values[index]-reference)
            # Float64 matrix product, prefix and interval sums. Same inputs,
            # independent reduction order; derived accumulation scale.
            allowance=16*120*120*np.finfo(float).eps*max(1.,reference,float(abs(Z).max()))/math.sqrt(float(np.diag(S).min()))
            if abs(difference)>allowance:raise ArithmeticError('Independent Gaussian maximum mismatch')
            discrepancies.append(dict(index=index,difference=difference,allowance=allowance))
        ranks.append(dict(condition=condition,rank=k,critical=r['critical'],calibration_failure=failure,
                          scalar_maximum_checks=discrepancies))
    coverage=[];eps=np.finfo(float).eps
    if len(summary['coverage'])!=12:raise AssertionError('Incomplete coverage family')
    for r in summary['coverage']:
        subset=frame[frame.phase.eq('coverage')&frame.method.eq(r['method'])&frame.condition.eq(r['condition'])&frame.interval.eq(r['interval'])]
        k=int(subset.covered.sum());n=len(subset)
        if (n,k)!=(r['n'],r['covered']) or n!=1000:raise AssertionError('Coverage count mismatch')
        tail=.05/(2*12)
        lower=0. if k==0 else optimize.brentq(lambda p:stats.binom.sf(k-1,n,p)-tail,0.,1.,xtol=eps,rtol=4*eps)
        upper=1. if k==n else optimize.brentq(lambda p:stats.binom.cdf(k,n,p)-tail,0.,1.,xtol=eps,rtol=4*eps)
        if max(abs(lower-r['simultaneous_cp_low']),abs(upper-r['simultaneous_cp_high']))>64*eps:
            raise ArithmeticError('Independent coverage bound mismatch')
        coverage.append(dict(method=r['method'],condition=r['condition'],interval=r['interval'],covered=k,n=n,
                             simultaneous_cp_low=lower,simultaneous_cp_high=upper))
    if len(summary['power'])!=48:raise AssertionError('Incomplete power reporting')
    for r in summary['power']:
        subset=frame[frame.phase.eq('power')&frame.method.eq(r['method'])&frame.condition.eq(r['condition'])&
                     frame.effect.eq(r['effect'])&frame.interval.eq(r['interval'])]
        if len(subset)!=200 or int(subset.covered.sum())!=r['covered'] or int(subset.excludes_zero.sum())!=r['zero_exclusions']:
            raise AssertionError('Power count mismatch')
    base=json.loads((a.task/'baseline/source_manifest.json').read_text());unchanged=0
    for path,h in base.items():
        if hashlib.sha256((a.task/'baseline'/path).read_bytes()).hexdigest()!=h:raise AssertionError('Baseline changed')
        if path.startswith(('python/','tests/','crates/')):
            if hashlib.sha256((a.task.parents[1]/path).read_bytes()).hexdigest()!=h:raise AssertionError('Core/test/native source changed')
            unchanged+=1
    report=dict(input_hashes_verified=len(manifest['input']['files']),original_intervals_preserved=10800,
                candidate_interval_selections_unchanged=10800,calibration=ranks,coverage=coverage,power_endpoints_recounted=48,
                baseline_files_preserved=len(base),core_test_native_unchanged=unchanged,
                scope='Known-covariance simulation and numerical audit; not biological calibration')
    (a.out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)


if __name__=='__main__':main()
