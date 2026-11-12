"""Conditional RNA-only convex stacking of the existing profile smoothers."""
import argparse
from fractions import Fraction
from functools import lru_cache
import importlib.util
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import linprog
from scipy.stats import pearsonr
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record
from haplodose.smoothing import predictive_rna_profile, _running_median


OLD = ('pseudobulk', 'pseudobulk_sm5', 'pseudobulk_mdl', 'predictive_smoothed', 'predictive_partition')
LOCKED = ('2313287_STOMACH', '42MGBA_CENTRAL_NERVOUS_SYSTEM')


@lru_cache(None)
def load(name):
    spec = importlib.util.spec_from_file_location('stacked_'+name, Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def convex_l1(P, target):
    """Simplex L1 fit with an independently evaluated feasible dual bound."""
    P, target = np.asarray(P, float), np.asarray(target, float)
    if P.ndim != 2 or target.shape != (len(P),) or not P.size or not np.isfinite(P).all() or not np.isfinite(target).all():
        raise ValueError('Nonempty finite aligned expert predictions and targets required')
    n, k = P.shape
    identity = sparse.eye(n, format='csr')
    A = sparse.vstack([sparse.hstack([sparse.csr_matrix(P), -identity]),
        sparse.hstack([-sparse.csr_matrix(P), -identity])], format='csr')
    eq = sparse.csr_matrix((np.ones(k), (np.zeros(k, int), np.arange(k))), shape=(1, k+n))
    result = linprog(np.r_[np.zeros(k), np.full(n, 1/n)], A_ub=A,
        b_ub=np.r_[target, -target], A_eq=eq, b_eq=[1.], bounds=(0., None), method='highs',
        options={'threads': 1})
    if not result.success:
        return None, dict(status='unresolved_solver', solver_status=int(result.status), message=result.message)
    positive = np.maximum(result.x[:k], 0.)
    total = math.fsum(map(float, positive))
    if not total > 0:
        raise ArithmeticError('No positive simplex mass')
    weights = positive/total
    prediction = P@weights
    primal = math.fsum(map(float, abs(prediction-target)))/n
    # min_w max_|s|<=1 s'(Pw-y)/n. Any clipped s gives a valid
    # lower bound min_j P_j's/n-y's/n, regardless of solver residuals.
    dual = np.clip(n*(result.ineqlin.marginals[n:]-result.ineqlin.marginals[:n]), -1., 1.)
    products = [math.fsum(float(a*b) for a,b in zip(P[:,j], dual)) for j in range(k)]
    lower = (min(products)-math.fsum(float(a*b) for a,b in zip(target, dual)))/n
    operations = 32*n*(k+1)+256
    gamma = operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)
    magnitude = max(1., float(np.max(abs(P))), float(np.max(abs(target))))
    rounding = gamma*magnitude
    if primal+rounding < lower-rounding:
        raise ArithmeticError('Primal-dual ordering contradicts derived rounding bounds')
    return weights, dict(status='estimated', weights=weights.tolist(), mean_absolute_error=primal,
        feasible_dual_lower=lower, optimization_gap_upper=max(0., primal-lower)+2*rounding,
        rounding_budget=rounding, simplex_sum=float(weights.sum()), solver_message=result.message,
        certificate='Feasible clipped L1 dual; not an inferential confidence interval')


def training_predictions(values, gene_bins, blocks, widths, folds, excluded):
    sizes = np.bincount(gene_bins, minlength=len(blocks))
    hold = (folds == excluded)&(sizes[gene_bins] >= 2)
    track = np.array([np.median(values[(~hold)&(gene_bins == b)]) for b in range(len(blocks))])
    medians = np.array([_running_median(track, blocks, w) for w in widths])
    segmented, decision = load('rna_predictive_partition').segmented(track, blocks)
    return medians, segmented, hold, decision


def select(values, gene_bins, blocks):
    values, gene_bins, blocks = np.asarray(values, float), np.asarray(gene_bins, int), np.asarray(blocks)
    baseline, record = predictive_rna_profile(values, gene_bins, blocks)
    if record['selected_width'] is None:
        return baseline, dict(status='insufficient_masked_gene_support', selected=False, baseline=record)
    widths = record['widths']
    folds = pd.Series(gene_bins).groupby(gene_bins, sort=False).cumcount().to_numpy()%10
    eligible = np.bincount(gene_bins, minlength=len(blocks))[gene_bins] >= 2
    P = np.full((len(values), len(widths)+1), np.nan)
    partition_ok, fold_records = True, []
    for f in range(10):
        medians, segmented, hold, decision = training_predictions(values, gene_bins, blocks, widths, folds, f)
        if not hold.any():continue
        P[np.ix_(hold, np.arange(len(widths)))] = medians[:, gene_bins[hold]].T
        if segmented is None:partition_ok = False
        else:P[hold, -1] = segmented[gene_bins[hold]]
        fold_records.append(dict(fold=f, targets=int(hold.sum()), partition=decision))
    full = np.array([np.median(values[gene_bins == b]) for b in range(len(blocks))])
    F = np.array([_running_median(full, blocks, w) for w in widths])
    segmented, full_partition = load('rna_predictive_partition').segmented(full, blocks)
    partition_ok = partition_ok and segmented is not None
    names = [f'median_width_{w}' for w in widths]
    if partition_ok:
        F = np.vstack([F, segmented]);names.append('piecewise_mean_mdl')
    else:P = P[:, :-1]
    if not np.isfinite(P[eligible]).all():raise ArithmeticError('Incomplete out-of-training prediction matrix')
    weights, fit = convex_l1(P[eligible], values[eligible])
    single_risk = float(np.min(np.mean(abs(P[eligible]-values[eligible,None]), axis=0)))
    selected = bool(weights is not None and fit['mean_absolute_error']+fit['optimization_gap_upper']+
        2*fit['rounding_budget'] < single_risk)
    if selected:
        output = weights@F
    else:
        # The installed median is the prespecified unresolved-decision fallback.
        output = baseline
    return output, dict(status=fit['status'], selected=selected, baseline=record, experts=names,
        best_single_expert_risk=single_risk, stacking=fit, eligible_genes=int(eligible.sum()),
        full_partition=full_partition, fold_partitions=fold_records,
        interpretation='Conditional descriptive RNA stack; no independent confirmation or count calibration')


def checks():
    P = np.array([[0.,2.],[0.,2.],[2.,4.]])
    y = np.array([1.,2.,3.])
    w, result = convex_l1(P,y)
    np.testing.assert_array_equal(w,[.5,.5])
    exact = Fraction(1,3)
    assert abs(result['mean_absolute_error']-float(exact)) <= result['rounding_budget']
    assert abs(result['feasible_dual_lower']-float(exact)) <= result['rounding_budget']
    # Exhaustive exact breakpoints independently minimize the two-expert L1.
    candidates = {Fraction(0),Fraction(1)}
    for a,b,t in zip(P[:,0],P[:,1],y):
        candidates.add(min(Fraction(1),max(Fraction(0),(Fraction(float(t))-Fraction(float(a)))/(Fraction(float(b))-Fraction(float(a))))))
    risks = [sum(abs((1-z)*Fraction(float(a))+z*Fraction(float(b))-Fraction(float(t))) for a,b,t in zip(P[:,0],P[:,1],y))/len(y) for z in candidates]
    assert min(risks) == exact
    permuted, alternate = convex_l1(P[::-1],y[::-1])
    np.testing.assert_array_equal(w,permuted)
    gbin = np.repeat(np.arange(7),4);blocks = np.array(['a']*4+['b']*3)
    values = np.repeat([0.,0.,1.,1.,10.,10.,11.],4)
    folds = pd.Series(gbin).groupby(gbin,sort=False).cumcount().to_numpy()%10
    before, segment, hold, _ = training_predictions(values,gbin,blocks,[1,3],folds,0)
    changed = values.copy();changed[hold] += 100.
    after, other, _, _ = training_predictions(changed,gbin,blocks,[1,3],folds,0)
    np.testing.assert_array_equal(before,after)
    if segment is not None:np.testing.assert_array_equal(segment,other)
    constant, decision = select(np.ones(len(values)),gbin,blocks)
    np.testing.assert_array_equal(constant,np.ones(7));assert not decision['selected']
    singleton, decision = select(np.array([1.,2.,3.]),np.arange(3),np.array(['a']*3))
    np.testing.assert_array_equal(singleton,[1.,2.,3.]);assert decision['status']=='insufficient_masked_gene_support'
    return dict(exact_breakpoint_optimum=str(exact),exact_half_weights=True,
        primal_dual_certificate=result,masked_base_predictions_invariant=True,
        constant_retains_original=True,singletons_abstain=True,row_permutation=True)


def panel(args):
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    prior = root/'dev/predictive_partition_20261109/kinker_v1'
    base = root/'dev/reference_weighted_pseudobulk_20261109/regression_v1'
    identities = json.loads((base/'predictions_locked.json').read_text())['cases']
    selected = identities if args.all_lines else [r for r in identities if r['case'] in LOCKED]
    assert len(selected) == (181 if args.all_lines else 2)
    if args.all_lines:
        assert 0 <= args.start < args.stop <= 181
        selected = selected[args.start:args.stop]
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        partition=Path(__file__).with_name('rna_predictive_partition.py'),
        partition_engine=Path(__file__).with_name('performance_candidates.py'),
        identities=base/'predictions_locked.json',old_scores=prior/'scores.tsv',
        truth=root/'results/kinker_input_v2_20261101_1551/depmap_cn.tsv.gz')
    if args.all_lines:inputs['regression_design'] = args.task/'REGRESSION_DESIGN.md'
    if args.all_lines:
        inputs['locked_gate'] = args.task/'locked_v2/report.json'
        assert json.loads(inputs['locked_gate'].read_text())['locked_two_case_gate']
    for r in selected:
        i = r['prediction_index'];inputs[f'RNA{i}'] = prior/f'predictions_case{i}.tsv.gz'
        inputs[f'bins{i}'] = base/f'bins_case{i}.tsv'
    run = start_run(args.out, inputs, dict(all_lines=args.all_lines, fold_count=10,
        experts='Original geometric median widths plus existing MDL', optimizer='simplex L1 with feasible dual'),
        dict(scope='Conditional prepared RNA prediction; former evaluation regression-only',
            biological_calibration=False,counts='CPM reconstruction; no inferential count model'),mode='rna_stacked_smoothing')
    try:
        numerical = checks();predictions = [];decisions = []
        for r in selected:
            i = r['prediction_index'];gene = pd.read_csv(inputs[f'RNA{i}'],sep='\t',float_precision='round_trip')
            bins = pd.read_csv(inputs[f'bins{i}'],sep='\t')
            assert gene.gene.is_unique and gene.case.eq(r['case']).all()
            value = gene.rna_log_ratio.to_numpy();gbin = gene.bin.to_numpy();blocks = bins.smoothing_block.to_numpy()
            unchanged, _ = predictive_rna_profile(value,gbin,blocks)
            np.testing.assert_array_equal(unchanged[gbin],gene.predictive_smoothed)
            fitted, decision = select(value,gbin,blocks)
            gene['stacked_smoothed'] = fitted[gbin]
            gene.to_csv(args.out/f'predictions_case{i}.tsv.gz',sep='\t',index=False)
            decisions.append(dict(case=r['case'],prediction_index=i,**decision));predictions.append(gene)
        write_json(dict(cases=decisions,files=[file_record(args.out/f"predictions_case{r['prediction_index']}.tsv.gz") for r in selected]),
            args.out/'RNA_DECISIONS_LOCKED.json')
        cn = pd.read_csv(inputs['truth'],sep='\t',index_col=0)
        assert cn.index.is_unique and cn.columns.is_unique
        rows = []
        methods = OLD+('stacked_smoothed',)
        for gene,r in zip(predictions,selected):
            assert len(r['query_DNA_line']) == 1
            weights = pd.Series(r['reference_DNA_weights'])
            g = gene[gene.gene.isin(cn.index)].copy()
            g['truth'] = cn.loc[g.gene,r['query_DNA_line'][0]].to_numpy()-cn.loc[g.gene,weights.index].to_numpy()@weights.to_numpy()
            assert np.isfinite(g[['truth']+list(methods)].to_numpy()).all()
            g['mb'] = g.chromosome+':'+((g.start.to_numpy()+g.end.to_numpy())//2//1_000_000).astype(str)
            mb = g.groupby('mb')[['truth']+list(methods)].mean()
            mb.to_csv(args.out/f"windows_case{r['prediction_index']}.tsv",sep='\t')
            for name in methods:rows.append(dict(case=r['case'],method=name,
                pearson=float(pearsonr(mb[name],mb.truth).statistic),windows=len(mb),common_genes=len(g)))
        scores = pd.DataFrame(rows);scores.to_csv(args.out/'scores.tsv',sep='\t',index=False)
        previous = pd.read_csv(inputs['old_scores'],sep='\t',float_precision='round_trip')
        old = scores[scores.method.isin(OLD)]
        before = previous[previous.case.isin([r['case'] for r in selected])]
        # Reparse the actual CSV to compare at the same frozen float values.
        after = pd.read_csv(args.out/'scores.tsv',sep='\t',float_precision='round_trip')
        pd.testing.assert_frame_equal(before.reset_index(drop=True),after[after.method.isin(OLD)].reset_index(drop=True),check_exact=True)
        matrix = scores.pivot(index='case',columns='method',values='pearson')
        paired = [dict(baseline=name,cases=len(matrix),mean=float((matrix.stacked_smoothed-matrix[name]).mean()),
            wins=int((matrix.stacked_smoothed>matrix[name]).sum()),losses=int((matrix.stacked_smoothed<matrix[name]).sum()),
            ties=int((matrix.stacked_smoothed==matrix[name]).sum())) for name in OLD]
        locked = matrix.loc[matrix.index.isin(LOCKED)]
        gate = len(locked) == 2 and all((locked.stacked_smoothed > locked[name]).all()
            for name in ('pseudobulk_sm5','predictive_smoothed','predictive_partition'))
        write_json(dict(cases=len(selected),checks=numerical,means=scores.groupby('method').pearson.mean().to_dict(),
            paired=paired,locked_two_case_gate=bool(gate),previous_scores_exact=len(old),
            experimental_stack_selected=sum(r['selected'] for r in decisions),independent_evaluation=False,
            biological_calibration=False,production_changed=False,scientific_criteria_changed=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU sparse HiGHS L1 stacking and matched-support DNA scoring')
    except BaseException as error:
        finish_run(args.out,run,error=error);raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--checks-only',action='store_true');parser.add_argument('--all-lines',action='store_true')
    parser.add_argument('--start',type=int,default=0);parser.add_argument('--stop',type=int,default=181)
    args = parser.parse_args();args.out.mkdir(exist_ok=False)
    if args.checks_only:write_json(checks(),args.out/'report.json')
    else:panel(args)


if __name__ == '__main__':main()
