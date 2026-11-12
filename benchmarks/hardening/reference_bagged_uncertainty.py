"""Conditional patient uncertainty for all fixed reference-bagging contrasts."""
import argparse
from fractions import Fraction
from functools import lru_cache
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, binom
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


@lru_cache(None)
def load(name):
    spec = importlib.util.spec_from_file_location('bag_uncertainty_'+name,
        Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    primary = load('reference_bagged_profile').PRIMARY
    inputs = dict(script=Path(__file__),design=args.task/'CONDITIONAL_DESIGN.md',
        scores=args.task/'score_v1/patient_scores.tsv',report=args.task/'score_v1/report.json',
        statistics=Path(__file__).with_name('statistics.py'),
        independent=Path(__file__).with_name('gene_factorial_precision_audit.py'))
    rec = start_run(args.out,inputs,dict(seed=2026111031,draws=20000,family=249),
        dict(scope='Approximate patient bootstrap CONDITIONAL on frozen reference ensemble/calibration; no full reference uncertainty'),
        mode='reference_bag_conditional_uncertainty')
    try:
        table = pd.read_csv(inputs['scores'],sep='\t',float_precision='round_trip')
        table = table[table.split.eq('evaluation')]
        methods = table.method.drop_duplicates().tolist()
        assert len(methods)==84 and primary in methods
        definitions,values = [],[]
        for name in methods:
            if name==primary:continue
            for endpoint in ('pearson','recall_loss','fpr_loss'):
                matrix = table.pivot(index='sample',columns='method',values=endpoint)
                assert len(matrix)==24
                direction = -1 if endpoint=='fpr_loss' else 1
                values.append((matrix[primary]-matrix[name]).to_numpy(float)*direction)
                definitions.append(dict(baseline=name,endpoint=endpoint,direction=direction))
        values = np.asarray(values)
        finite = np.isfinite(values)
        assert values.shape==(249,24)
        count = finite.sum(axis=1)
        point = np.array([math.fsum(map(float,row[np.isfinite(row)]))/int(n)
                          for row,n in zip(values,count)])
        rng = np.random.default_rng(2026111031)
        indices = rng.integers(0,24,size=(20000,24))
        weights = np.array([np.bincount(row,minlength=24)for row in indices])
        denominator = weights@finite.astype(int).T
        numerators = weights@np.where(finite,values,0).T
        draws = np.full(numerators.shape,np.nan)
        np.divide(numerators,denominator,out=draws,where=denominator>0)
        # Independent scalar sums check numerical matrix accumulation on actual
        # draws, with an operation/magnitude rounding allowance.
        eps = np.finfo(float).eps
        gamma = 128*24*eps/(1-128*24*eps)
        error = 0.
        for i in (0,1,19999):
            for j,row in enumerate(values):
                selected = row[indices[i]]
                selected = selected[np.isfinite(selected)]
                if len(selected):
                    exact = math.fsum(map(float,selected))/len(selected)
                    discrepancy = abs(exact-draws[i,j])
                    assert discrepancy<=gamma*max(1.,float(np.max(abs(selected))))
                    error = max(error,discrepancy)
        np.savez_compressed(args.out/'draws.npz',indices=indices,draws=draws,paired=denominator)
        reference = load('gene_factorial_precision_audit')
        family = len(definitions)
        mc_tail = .05/(2*3*family)
        rows = []
        for j,definition in enumerate(definitions):
            b = draws[:,j]
            b = b[np.isfinite(b)]
            delta,n = point[j],len(b)
            testable = bool(n==20000 and count[j]>=2 and np.ptp(b)>0)
            row = dict(**definition,paired_patients=int(count[j]),effect=float(delta),
                valid_draws=n,testable=testable,p=1.,p_mc_lower=0.,p_mc_upper=1.,
                family_ci_low=None,family_ci_high=None,
                family_ci_low_mc=[-np.inf,np.inf],family_ci_high_mc=[-np.inf,np.inf])
            if n:
                errors = b-delta
                ordered = np.sort(errors)
                for key,q in (('family_ci_low',1-.05/(2*family)),('family_ci_high',.05/(2*family))):
                    value = float(delta-np.quantile(errors,q))
                    independent,guard = reference.basic_endpoint(ordered,delta,q)
                    assert abs(value-independent)<=guard
                    row[key] = value
                    lo = int(binom.ppf(mc_tail,n,q))
                    hi = int(binom.isf(mc_tail,n,q))+1
                    assert binom.cdf(lo-1,n,q)<=mc_tail and binom.sf(hi-1,n,q)<=mc_tail
                    lower = -np.inf if lo==0 else float(ordered[lo-1])
                    upper = np.inf if hi>n else float(ordered[hi-1])
                    row[key+'_mc'] = [float(delta-upper),float(delta-lower)]
                if testable:
                    k = int(np.count_nonzero(abs(errors)>=abs(delta)))
                    row['p'] = float(Fraction(k+1,n+1))
                    row['p_mc_lower'] = 0. if k==0 else float(beta.ppf(mc_tail,k,n-k+1))
                    row['p_mc_upper'] = 1. if k==n else float(beta.isf(mc_tail,k+1,n-k))
            rows.append(row)
        stats = load('statistics')
        for key in ('p','p_mc_lower','p_mc_upper'):
            adjusted = stats.holm([r[key]for r in rows])
            independent = reference.holm([r[key]for r in rows])
            np.testing.assert_array_equal(adjusted,independent)
            for row,value in zip(rows,adjusted):row[key+'_holm'] = float(value)
        for row in rows:
            row['conditional_gate_passed'] = bool(row['testable']and row['effect']>0
                and row['p_mc_upper_holm']<.05 and row['family_ci_low_mc'][0]>0)
            row['numerically_unresolved'] = bool(row['testable']and(
                row['p_mc_lower_holm']<=.05<=row['p_mc_upper_holm']or
                row['family_ci_low_mc'][0]<=0<=row['family_ci_low_mc'][1]))
        pearson = [r for r in rows if r['endpoint']=='pearson']
        write_json(dict(family=249,draws=20000,comparisons=rows,
            scalar_draw_check_maximum_error=error,independent_Holm_exact=True,
            profile_conditional_gate=all(r['conditional_gate_passed']for r in pearson),
            joint_conditional_gate=all(r['conditional_gate_passed']for r in rows),
            profile_comparisons_passed=sum(r['conditional_gate_passed']for r in pearson),
            all_comparisons_passed=sum(r['conditional_gate_passed']for r in rows),
            reference_uncertainty_included=False,calibration_uncertainty_included=False,
            coverage_validated=False,method_promoted=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU paired patient bootstrap,249 family,exact-binomial Monte Carlo brackets')
    except BaseException as error:
        finish_run(args.out,rec,error=error)
        raise


if __name__=='__main__':main()
