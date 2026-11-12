"""Research physical-window weighting of an otherwise unchanged RNA stack."""
import argparse
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from haplodose import gene_segments
from haplodose.stacking import _convex_l1
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


LOCKED = ('2313287_STOMACH', '42MGBA_CENTRAL_NERVOUS_SYSTEM')
OLD = ('pseudobulk', 'pseudobulk_sm5', 'pseudobulk_mdl', 'predictive_smoothed',
       'predictive_partition', 'stacked_smoothed', 'gene_laplace_alone', 'gene_laplace_stack')


def factors(windows):
    _, index, counts = np.unique(windows, return_inverse=True, return_counts=True)
    if not len(counts):
        raise ValueError('Identified target windows required')
    # _convex_l1 averages rows. This produces equal total loss weight per
    # represented physical window, then equal gene weight within that window.
    return len(index)/(len(counts)*counts[index])


def numerical_checks():
    P = np.tile([0., 4.], (10, 1))
    y = np.r_[np.zeros(8), 2., 2.]
    windows = np.r_[np.zeros(8, int), 1, 2]
    a = factors(windows)
    old, old_record = _convex_l1(P, y)
    new, record = _convex_l1(P*a[:, None], y*a)
    assert old is not None and new is not None
    np.testing.assert_array_equal(old, [1., 0.])
    # Independent exact minimization: an L1 two-expert optimum occurs at an
    # endpoint or a residual breakpoint. Enumerate with rational coefficients.
    candidates = {Fraction(0), Fraction(1)}
    for row, target in zip(P, y):
        if row[0] != row[1]:
            t = (Fraction(float(target))-Fraction(float(row[1])))/(
                 Fraction(float(row[0]))-Fraction(float(row[1])))
            if 0 <= t <= 1:
                candidates.add(t)
    count = {w: int(np.sum(windows == w)) for w in set(windows)}
    def exact_loss(t):
        return sum(abs(t*Fraction(float(row[0]))+(1-t)*Fraction(float(row[1]))
                       -Fraction(float(target)))/count[w]
                   for row, target, w in zip(P, y, windows))/len(count)
    minimum = min(exact_loss(t) for t in candidates)
    assert minimum == Fraction(2, 3)
    error = abs(record['mean_absolute_error']-float(minimum))
    assert error <= record['optimization_gap_upper']+record['rounding_budget']
    assert abs(new[0]-.5) <= record['rounding_budget']
    # Duplicating every gene in the dense window preserves the target loss,
    # unlike giving each observed gene a separate equally weighted target.
    indices = np.r_[np.arange(8), np.arange(10)]
    weight = factors(windows[indices])
    duplicated, duplicate_record = _convex_l1(P[indices]*weight[:, None], y[indices]*weight)
    assert duplicated is not None
    assert np.max(abs(duplicated-new)) <= (
        record['rounding_budget']+duplicate_record['rounding_budget'])
    uniform = factors(np.arange(len(y)))
    np.testing.assert_array_equal(uniform, np.ones(len(y)))
    same, _ = _convex_l1(P*uniform[:, None], y*uniform)
    np.testing.assert_array_equal(same, old)
    return dict(exact_two_expert_minimum='2/3', exact_breakpoint_error=error,
        unique_balanced_first_weight=.5, unweighted_first_weight=float(old[0]),
        dense_window_replication_invariant=True, uniform_reduces_to_original=True,
        numerical_budget=record['optimization_gap_upper']+record['rounding_budget'],
        scientific_thresholds_changed=False)


def fit(values, gene_bins, blocks, windows):
    if sys.getprofile() is not None:
        raise RuntimeError('Expert capture requires an otherwise free observation hook')
    captured = {}
    def observe(frame, event, result):
        if event == 'return' and frame.f_code is gene_segments._select.__code__:
            local = frame.f_locals
            if 'P' in local and 'eligible' in local and 'F' in local:
                captured.update(P=local['P'].copy(), eligible=local['eligible'].copy(),
                    bank=np.vstack([local['F'], local['full']]),
                    widths=local['widths'], previous=local['previous'])
    sys.setprofile(observe)
    try:
        baseline, decision, _ = gene_segments._select(values, gene_bins, blocks)
    finally:
        sys.setprofile(None)
    if not captured:
        return baseline, baseline, dict(selected=False, reason='existing_expert_support_unavailable'), {}
    eligible = captured['eligible']
    P, y = captured['P'][eligible], np.asarray(values)[eligible]
    a = factors(np.asarray(windows)[eligible])
    if decision['selected']:
        original_weights = np.asarray(decision['stacking']['weights'])
    elif captured['previous']['selected']:
        original_weights = np.r_[captured['previous']['stacking']['weights'], 0.]
    else:
        original_weights = np.zeros(P.shape[1])
        index = captured['widths'].index(captured['previous']['baseline']['selected_width'])
        original_weights[index] = 1.
    np.testing.assert_array_equal(original_weights@captured['bank'], baseline)
    old_risk = math.fsum(map(float, a*abs(P@original_weights-y)))/len(y)
    weights, record = _convex_l1(P*a[:, None], y*a)
    operations = 32*len(y)*(P.shape[1]+4)+256
    eps = np.finfo(float).eps
    if operations*eps >= 1:
        raise ArithmeticError('Weighted risk arithmetic cannot be resolved')
    gamma = operations*eps/(1-operations*eps)
    budget = gamma*max(1., float(np.max(abs(P))), float(np.max(abs(y))))*max(1., float(a.max()))
    chosen = bool(weights is not None and record['mean_absolute_error']+
        record['optimization_gap_upper']+2*record['rounding_budget'] < old_risk-budget)
    output = weights@captured['bank'] if chosen else baseline
    diagnostics = dict(selected=chosen, original_balanced_risk=old_risk,
        original_risk_rounding=budget, weighted_fit=record,
        targets=len(y), physical_windows=len(set(np.asarray(windows)[eligible])),
        minimum_target_factor=float(a.min()), maximum_target_factor=float(a.max()),
        original_weights=original_weights.tolist(), experts=decision['experts'],
        units='natural_log_relative_RNA_rate', calibration_status='unvalidated',
        estimand='Mean within-window gene absolute prediction error,equally weighted over represented1Mb windows',
        independent_confirmation=False)
    arrays = dict(predictions=P, target=y, target_factor=a,
        target_windows=np.asarray(windows)[eligible], bank=captured['bank'],
        original_weights=original_weights, candidate_weights=weights if weights is not None else [],
        baseline=baseline, candidate=output)
    return output, baseline, diagnostics, arrays


def panel(args):
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    prior = root/'dev/gene_laplace_profile_20261110'
    base = root/'dev/reference_weighted_pseudobulk_20261109/regression_v1'
    identities = json.loads((base/'predictions_locked.json').read_text())['cases']
    selected = [r for r in identities if r['case'] in LOCKED]
    assert tuple(r['case'] for r in selected) == LOCKED
    numerics = args.task/'numerics_v1/report.json'
    assert json.loads(numerics.read_text())['dense_window_replication_invariant']
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md', numerics=numerics,
        identities=base/'predictions_locked.json', previous_scores=prior/'sc_v1/scores.tsv',
        truth=root/'results/kinker_input_v2_20261101_1551/depmap_cn.tsv.gz')
    for row in selected:
        i = row['prediction_index']
        inputs[f'RNA{i}'] = prior/f'sc_v1/predictions_case{i}.tsv.gz'
        inputs[f'bins{i}'] = base/f'bins_case{i}.tsv'
    run = start_run(args.out, inputs, dict(window_bp=1000000, expert_bank_unchanged=True),
        dict(scope='Two locked exposed profile cases;not independent evaluation or count calibration'),
        mode='window_balanced_stack')
    try:
        predictions, decisions = [], []
        for row in selected:
            i = row['prediction_index']
            gene = pd.read_csv(inputs[f'RNA{i}'], sep='\t', float_precision='round_trip')
            bins = pd.read_csv(inputs[f'bins{i}'], sep='\t')
            assert gene.gene.is_unique and gene.case.eq(row['case']).all()
            windows = gene.chromosome+':'+((gene.start+gene.end)//2//1000000).astype(str)
            candidate, baseline, decision, arrays = fit(gene.rna_log_ratio.to_numpy(),
                gene.bin.to_numpy(), bins.smoothing_block.to_numpy(), windows.to_numpy())
            np.testing.assert_array_equal(baseline[gene.bin.to_numpy()], gene.gene_laplace_stack)
            gene['window_balanced_stack'] = candidate[gene.bin.to_numpy()]
            file = args.out/f'predictions_case{i}.tsv.gz'
            gene.to_csv(file, sep='\t', index=False)
            np.savez_compressed(args.out/f'fit_case{i}.npz', **arrays)
            decisions.append(dict(case=row['case'], prediction_index=i, **decision))
            predictions.append(gene)
        write_json(dict(cases=decisions, outputs=[file_record(args.out/f"predictions_case{r['prediction_index']}.tsv.gz")
            for r in selected], DNA_parsed=False), args.out/'RNA_LOCKED.json')
        cn = pd.read_csv(inputs['truth'], sep='\t', index_col=0)
        assert cn.index.is_unique and cn.columns.is_unique
        rows = []
        for row, gene in zip(selected, predictions):
            ref = pd.Series(row['reference_DNA_weights'])
            g = gene[gene.gene.isin(cn.index)].copy()
            g['truth'] = cn.loc[g.gene, row['query_DNA_line'][0]].to_numpy()-cn.loc[g.gene, ref.index].to_numpy()@ref.to_numpy()
            g['mb'] = g.chromosome+':'+((g.start.to_numpy()+g.end.to_numpy())//2//1000000).astype(str)
            mb = g.groupby('mb')[['truth']+list(OLD)+['window_balanced_stack']].mean()
            for method in OLD+('window_balanced_stack',):
                finite = np.isfinite(mb[[method, 'truth']]).all().all() and mb[method].nunique() > 1
                score = float(pearsonr(mb[method], mb.truth).statistic) if finite else np.nan
                rows.append(dict(case=row['case'], method=method, pearson=score,
                                 windows=len(mb), common_genes=len(g)))
        scores = pd.DataFrame(rows)
        previous = pd.read_csv(inputs['previous_scores'], sep='\t', float_precision='round_trip')
        pd.testing.assert_frame_equal(previous, scores[scores.method.isin(OLD)].reset_index(drop=True), check_exact=True)
        scores.to_csv(args.out/'scores.tsv', sep='\t', index=False)
        matrix = scores.pivot(index='case', columns='method', values='pearson')
        required = ('gene_laplace_stack', 'pseudobulk_sm5', 'pseudobulk_mdl')
        gate = bool(all((matrix.window_balanced_stack > matrix[name]).all() for name in required))
        write_json(dict(gate=gate, cases=2, paired=matrix.reset_index().to_dict('records'),
            RNA_decisions=decisions, previous_scores_exact=len(previous), expert_bank_unchanged=True,
            independent_evaluation=False, count_calibration=False, production_changed=False,
            decision='Permit separately frozen full regression' if gate else
                     'Reject expansion;no window-size or target-weight tuning'), args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU same expert bank/simplex L1 with physical-window loss weights')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--mode', choices=('numerics', 'panel'), required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    if args.mode == 'numerics':
        write_json(numerical_checks(), args.out/'report.json')
    else:
        panel(args)


if __name__ == '__main__':
    main()
