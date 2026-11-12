"""Research averaging of the existing gene-Laplace segmentation objective."""
import argparse
from fractions import Fraction
from itertools import product
import json
import math
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import pearsonr
from haplodose import gene_segments
from haplodose.stacking import _convex_l1
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


LOCKED = ('2313287_STOMACH', '42MGBA_CENTRAL_NERVOUS_SYSTEM')
OLD = ('pseudobulk', 'pseudobulk_sm5', 'pseudobulk_mdl', 'predictive_smoothed',
       'predictive_partition', 'stacked_smoothed', 'gene_laplace_alone', 'gene_laplace_stack')


def average(costs, medians, unit, denominator, penalty):
    """Sum over every partition; no thresholding, truncation or path sampling."""
    n = len(costs)
    logw, effects = np.full((n, n), -np.inf), np.zeros((n, n))
    for lo in range(n):
        for hi in range(lo, n):
            logw[lo, hi] = -.5*float(Fraction(costs[lo][hi]+penalty, unit))
            effects[lo, hi] = float(Fraction(medians[lo][hi], 2*denominator))
    if not np.isfinite(logw[np.triu_indices(n)]).all():
        raise ArithmeticError('Finite interval scores required')
    forward, backward = np.full(n+1, -np.inf), np.full(n+1, -np.inf)
    forward[0], backward[-1] = 0., 0.
    for hi in range(1, n+1):
        forward[hi] = logsumexp(forward[:hi]+logw[:hi, hi-1])
    for lo in range(n-1, -1, -1):
        backward[lo] = logsumexp(logw[lo, lo:]+backward[lo+1:])
    changes, mass_changes = [[] for _ in range(n+1)], [[] for _ in range(n+1)]
    for lo in range(n):
        for hi in range(lo+1, n+1):
            mass = math.exp(forward[lo]+logw[lo, hi-1]+backward[hi]-forward[-1])
            value = mass*effects[lo, hi-1]
            changes[lo].append(value); changes[hi].append(-value)
            mass_changes[lo].append(mass); mass_changes[hi].append(-mass)
    delta = [math.fsum(v) for v in changes]
    mass_delta = [math.fsum(v) for v in mass_changes]
    fitted = np.array([math.fsum(delta[:hi]) for hi in range(1, n+1)])
    coverage = np.array([math.fsum(mass_delta[:hi]) for hi in range(1, n+1)])
    operations = 128*n*n+256
    eps = np.finfo(float).eps
    assert operations*eps < 1
    gamma = operations*eps/(1-operations*eps)
    log_scale = max(1., float(np.max(abs(logw[np.triu_indices(n)]))),
                    float(np.max(abs(forward))), float(np.max(abs(backward))))
    mass_budget = gamma*log_scale
    effect_budget = mass_budget*max(1., float(np.max(abs(effects))))
    assert np.max(abs(coverage-1)) <= mass_budget
    return fitted, dict(log_normalizer=float(forward[-1]),
        coverage_maximum_error=float(np.max(abs(coverage-1))),
        probability_arithmetic_budget=mass_budget, effect_arithmetic_budget=effect_budget,
        bins=n, intervals=n*(n+1)//2,
        interpretation='Normalized Gibbs weights on fitted segment losses;not calibrated boundary or DNA probabilities')


def fit(values, gene_bins, blocks):
    if sys.getprofile() is not None:
        raise RuntimeError('Read-only interval capture requires a free profile hook')
    pending, captured, folds = {}, {}, {}
    # Reuse the actual exact integer interval objectives in the old kernel.
    # Version-pinned return observation changes neither kernel nor its result.
    def observe(frame, event, result):
        if event != 'return':
            return
        loc, caller = frame.f_locals, frame.f_back
        if frame.f_code is gene_segments._block_partition.__code__:
            assert caller.f_code is gene_segments._partition_profile.__code__
            prediction, record = average(loc['costs'], loc['medians'], loc['unit'],
                                          loc['denominator'], loc['penalty_units'])
            pending.setdefault(id(caller), []).append((caller.f_locals['lo'],
                caller.f_locals['hi'], prediction, record))
        elif frame.f_code is gene_segments._partition_profile.__code__:
            parts = pending.pop(id(frame), [])
            if result is None or result[0] is None:
                assert not parts
                return
            assert caller.f_code is gene_segments._select.__code__ and parts
            output = np.full(len(result[0]), np.nan)
            for lo, hi, prediction, _ in parts:
                assert np.isnan(output[lo:hi]).all()
                output[lo:hi] = prediction
            assert np.isfinite(output).all()
            key = caller.f_locals.get('fold', 'full')
            assert key not in folds
            folds[key] = dict(profile=output, records=[p[3] for p in parts])
        elif frame.f_code is gene_segments._select.__code__:
            if all(key in loc for key in ('P', 'eligible', 'F', 'full')):
                captured.update(P=loc['P'].copy(), eligible=loc['eligible'].copy(),
                    bank=np.vstack([loc['F'], loc['full']]), previous=loc['previous'],
                    widths=loc['widths'])
    sys.setprofile(observe)
    try:
        baseline, old_decision, _ = gene_segments._select(values, gene_bins, blocks)
    finally:
        sys.setprofile(None)
    assert not pending
    if not captured:
        return baseline, baseline, dict(selected=False, reason='old_expert_support_unavailable'), {}
    values, gene_bins = np.asarray(values), np.asarray(gene_bins)
    ranks = pd.Series(gene_bins).groupby(gene_bins, sort=False).cumcount().to_numpy()%10
    eligible, old_P = captured['eligible'], captured['P']
    column = np.full(len(values), np.nan)
    for fold in sorted(key for key in folds if key != 'full'):
        use = eligible & (ranks == fold)
        column[use] = folds[fold]['profile'][gene_bins[use]]
    assert 'full' in folds and np.isfinite(column[eligible]).all()
    previous = captured['previous']
    if old_decision['selected']:
        weights = np.asarray(old_decision['stacking']['weights'])
        old_prediction = old_P[eligible]@weights
        np.testing.assert_array_equal(weights@captured['bank'], baseline)
    elif previous['selected']:
        weights = np.asarray(previous['stacking']['weights'])
        old_prediction = old_P[eligible, :len(weights)]@weights
        np.testing.assert_array_equal(weights@captured['bank'][:len(weights)], baseline)
    else:
        index = captured['widths'].index(previous['baseline']['selected_width'])
        old_prediction = old_P[eligible, index]
        np.testing.assert_array_equal(captured['bank'][index], baseline)
    target = values[eligible]
    old_risk = math.fsum(map(float, abs(old_prediction-target)))/len(target)
    P = np.c_[old_P[eligible], column[eligible]]
    weights, record = _convex_l1(P, target)
    averaging_budget = max(rec['effect_arithmetic_budget'] for key, value in folds.items()
                           for rec in value['records'])
    unit = (32*len(target)*(P.shape[1]+4)+256)*np.finfo(float).eps
    assert unit < 1
    old_risk_budget = unit/(1-unit)*max(1., float(np.max(abs(P))), float(np.max(abs(target))))
    chosen = bool(weights is not None and record['mean_absolute_error']+
        record['optimization_gap_upper']+2*record['rounding_budget']+averaging_budget
        < old_risk-old_risk_budget)
    standalone = folds['full']['profile']
    bank = np.vstack([captured['bank'], standalone])
    candidate = weights@bank if chosen else baseline
    full_budget = max(rec['effect_arithmetic_budget'] for rec in folds['full']['records'])
    dot_unit = (4*len(bank)+32)*np.finfo(float).eps
    dot_error = dot_unit/(1-dot_unit)*max(1., float(np.max(abs(bank))))
    diagnostics = dict(selected=chosen, new_expert_weight=float(weights[-1]) if weights is not None else None,
        old_risk=old_risk, old_risk_rounding=old_risk_budget, averaging_rounding=averaging_budget,
        standalone_profile_rounding=full_budget,
        candidate_profile_rounding=float(weights[-1])*full_budget+dot_error if chosen else 0.,
        profile_rounding_scope='Conditional on recorded float stacking weights;no uniqueness or sampling interval',
        stacking=record, fold_records={str(k):v['records'] for k,v in folds.items()},
        expert='All-partition average under existing gene-Laplace scores',
        calibration_status='unvalidated', units='natural_log_relative_RNA_rate')
    arrays = dict(old_predictions=old_P[eligible], average_predictions=column[eligible],
        target=target, bank=bank, baseline=baseline, candidate=candidate, standalone=standalone)
    return candidate, standalone, diagnostics, arrays


def checks():
    # With zero scores every partition has exactly equal rational mass.
    n = 4
    costs = [[0]*n for _ in range(n)]
    medians = [[2*(lo+hi) for hi in range(n)] for lo in range(n)]
    actual, record = average(costs, medians, 1, 1, 0)
    reference = [Fraction(0)]*n
    for cuts in product((False, True), repeat=n-1):
        edges = [0]+[i+1 for i,v in enumerate(cuts) if v]+[n]
        for lo, hi in zip(edges[:-1], edges[1:]):
            for j in range(lo, hi):
                reference[j] += Fraction(medians[lo][hi-1], 2*2**(n-1))
    discrepancy = float(np.max(abs(actual-np.array(list(map(float, reference))))))
    assert discrepancy <= record['effect_arithmetic_budget']
    # Capture actual median/loss tables and compare the entire fitted average
    # against independent exhaustive partition construction from raw genes.
    groups = [np.array([0., 0.]), np.array([1., 7.]), np.array([2., 3.]), np.array([9., 11.])]
    def block(groups, scale):
        captured = {}
        def observe(frame, event, result):
            if event == 'return' and frame.f_code is gene_segments._block_partition.__code__:
                v = frame.f_locals
                captured['result'] = average(v['costs'], v['medians'], v['unit'],
                    v['denominator'], v['penalty_units'])
        assert sys.getprofile() is None
        sys.setprofile(observe)
        try:
            gene_segments._block_partition(groups, scale)
        finally:
            sys.setprofile(None)
        return captured['result']
    scale = .75/math.log(2.)
    fitted, error = block(groups, scale)
    paths = []
    penalty = Fraction(3*math.log(sum(map(len, groups))))
    for cuts in product((False, True), repeat=3):
        edges = [0]+[i+1 for i,v in enumerate(cuts) if v]+[4]
        objective, prediction = penalty*(len(edges)-1), [None]*4
        for lo, hi in zip(edges[:-1], edges[1:]):
            observed = sorted(Fraction(float(v)) for g in groups[lo:hi] for v in g)
            median = (observed[(len(observed)-1)//2]+observed[len(observed)//2])/2
            objective += 2*sum(abs(v-median) for v in observed)/Fraction(scale)
            prediction[lo:hi] = [float(median)]*(hi-lo)
        paths.append((math.exp(-.5*float(objective)), prediction))
    total = math.fsum(w for w,_ in paths)
    reference = np.array([math.fsum(w*p[j] for w,p in paths)/total for j in range(4)])
    enumeration_error = float(np.max(abs(fitted-reference)))
    assert enumeration_error <= error['effect_arithmetic_budget']
    shifted, shift_error = block([g+2**20 for g in groups], scale)
    assert np.max(abs(shifted-(fitted+2**20))) <= shift_error['effect_arithmetic_budget']+error['effect_arithmetic_budget']
    scaled, scale_error = block([g*8 for g in groups], scale*8)
    assert np.max(abs(scaled-fitted*8)) <= scale_error['effect_arithmetic_budget']+8*error['effect_arithmetic_budget']
    singleton, _ = average([[0]], [[6]], 1, 1, 0)
    np.testing.assert_array_equal(singleton, [3.])
    values = np.tile([-.2,.1,.3,-.1,.2,-.3,.4,-.4,.5,-.5], 4)+np.repeat([0.,.5,0.,-.5],10)
    bins, blocks = np.repeat(np.arange(4), 10), np.zeros(4, int)
    _, _, _, first = fit(values, bins, blocks)
    changed = values.copy(); changed[0] += 64
    _, _, _, second = fit(changed, bins, blocks)
    assert first['average_predictions'][0] == second['average_predictions'][0]
    return dict(independent_enumeration_passed=True, rational_uniform_error=discrepancy,
        exhaustive_gene_loss_error=enumeration_error, effect_rounding_budget=error['effect_arithmetic_budget'],
        shift_and_scale_invariance=True, singleton=True, target_exclusion_exact=True)


def panel(args):
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    prior = root/'dev/gene_laplace_profile_20261110'
    base = root/'dev/reference_weighted_pseudobulk_20261109/regression_v1'
    selected = [r for r in json.loads((base/'predictions_locked.json').read_text())['cases'] if r['case'] in LOCKED]
    assert tuple(r['case'] for r in selected) == LOCKED
    numerics = args.task/'numerics_v1/report.json'
    assert json.loads(numerics.read_text())['independent_enumeration_passed']
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md', numerics=numerics,
        identities=base/'predictions_locked.json', previous_scores=prior/'sc_v1/scores.tsv',
        old_kernel=Path(gene_segments.__file__),
        truth=root/'results/kinker_input_v2_20261101_1551/depmap_cn.tsv.gz')
    for row in selected:
        i = row['prediction_index']
        inputs[f'RNA{i}'], inputs[f'bins{i}'] = prior/f'sc_v1/predictions_case{i}.tsv.gz', base/f'bins_case{i}.tsv'
    run = start_run(args.out, inputs, dict(cases=LOCKED, new_experts=1, old_loss_unchanged=True),
        dict(scope='Two exposed profile mechanisms;no independent evaluation or calibration'),
        mode='partition_averaged_stack')
    try:
        predictions, decisions = [], []
        for row in selected:
            i = row['prediction_index']
            genes = pd.read_csv(inputs[f'RNA{i}'], sep='\t', float_precision='round_trip')
            bins = pd.read_csv(inputs[f'bins{i}'], sep='\t')
            assert genes.gene.is_unique and genes.case.eq(row['case']).all()
            candidate, standalone, decision, arrays = fit(genes.rna_log_ratio.to_numpy(),
                genes.bin.to_numpy(), bins.smoothing_block.to_numpy())
            assert arrays
            np.testing.assert_array_equal(arrays['baseline'][genes.bin.to_numpy()], genes.gene_laplace_stack)
            genes['partition_average_stack'] = candidate[genes.bin.to_numpy()]
            genes['partition_average_alone'] = standalone[genes.bin.to_numpy()]
            genes.to_csv(args.out/f'predictions_case{i}.tsv.gz', sep='\t', index=False)
            np.savez_compressed(args.out/f'fit_case{i}.npz', **arrays)
            decisions.append(dict(case=row['case'], **decision))
            predictions.append(genes)
        write_json(dict(cases=decisions, outputs=[file_record(args.out/f"predictions_case{r['prediction_index']}.tsv.gz")
            for r in selected], DNA_parsed=False), args.out/'RNA_LOCKED.json')
        cn = pd.read_csv(inputs['truth'], sep='\t', index_col=0)
        assert cn.index.is_unique and cn.columns.is_unique
        rows, score_numerics = [], {}
        new = ('partition_average_alone', 'partition_average_stack')
        for row, genes in zip(selected, predictions):
            ref = pd.Series(row['reference_DNA_weights'])
            g = genes[genes.gene.isin(cn.index)].copy()
            g['truth'] = cn.loc[g.gene,row['query_DNA_line'][0]].to_numpy()-cn.loc[g.gene,ref.index].to_numpy()@ref.to_numpy()
            g['mb'] = g.chromosome+':'+((g.start.to_numpy()+g.end.to_numpy())//2//1000000).astype(str)
            mb = g.groupby('mb')[['truth']+list(OLD)+list(new)].mean()
            for method in OLD+new:
                finite = np.isfinite(mb[[method,'truth']]).all().all() and mb[method].nunique() > 1
                rows.append(dict(case=row['case'], method=method,
                    pearson=float(pearsonr(mb[method],mb.truth).statistic) if finite else np.nan,
                    windows=len(mb),common_genes=len(g)))
            decision = next(d for d in decisions if d['case'] == row['case'])
            centred = mb.partition_average_stack.to_numpy()-mb.partition_average_stack.mean()
            norm = float(np.linalg.norm(centred))
            perturbation = math.sqrt(len(mb))*decision['candidate_profile_rounding']
            # Centering is a contraction. Normalize y+e versus y, with
            # ||e||<=sqrt(windows)*max_bin_error; apply Cauchy-Schwarz to r.
            bound = 2*perturbation/(norm-perturbation) if norm > perturbation else np.inf
            score_numerics[row['case']] = dict(conditional_projection_error_bound=bound,
                centred_prediction_norm=norm, perturbation_norm_bound=perturbation,
                scope='Only export arithmetic at fixed fitted weights;not biological or weight-estimation uncertainty')
        scores = pd.DataFrame(rows)
        old = pd.read_csv(inputs['previous_scores'], sep='\t', float_precision='round_trip')
        pd.testing.assert_frame_equal(old,scores[scores.method.isin(OLD)].reset_index(drop=True),check_exact=True)
        scores.to_csv(args.out/'scores.tsv',sep='\t',index=False)
        matrix = scores.pivot(index='case',columns='method',values='pearson')
        required = ('gene_laplace_stack','pseudobulk_sm5','pseudobulk_mdl')
        margins = pd.Series({name:v['conditional_projection_error_bound'] for name,v in score_numerics.items()})
        gate = bool(all(((matrix.partition_average_stack-matrix[m]) > margins).all() for m in required)
            and all(d['selected'] and d['new_expert_weight'] > 0 for d in decisions))
        write_json(dict(gate=gate,cases=2,paired=matrix.reset_index().to_dict('records'),
            RNA_decisions=decisions,previous_scores_exact=len(old),independent_evaluation=False,
            score_numerics=score_numerics,
            count_calibration=False,production_changed=False,
            decision='Permit separately frozen full regression' if gate else
                'Reject expansion;no penalty,temperature,fold or case tuning'),args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(args.out,run,backend='CPU exact old interval losses plus all-partition Gibbs averaging')
    except BaseException as error:
        finish_run(args.out,run,error=error)
        raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--mode',choices=('numerics','panel'),required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    if args.mode=='numerics':
        write_json(checks(),args.out/'report.json')
    else:
        panel(args)


if __name__=='__main__':
    main()
