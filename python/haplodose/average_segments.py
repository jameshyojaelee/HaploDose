"""Descriptive averaging of gene-Laplace partitions in the RNA profile stack.

Gibbs weights on fitted working losses are not calibrated boundary probabilities.
Local callbacks reuse exact interval tables without global execution hooks.
"""
from fractions import Fraction
import math

import numpy as np
import pandas as pd
from scipy.special import logsumexp

from . import gene_segments
from .stacking import _convex_l1


def _average(costs, medians, unit, denominator, penalty):
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


def _fit(values, gene_bins, blocks, *, solver='primal'):
    """Return the conditional candidate, ablation, diagnostics and fit arrays."""
    folds = {}
    def partition(key):
        if key in folds:
            raise ArithmeticError('Duplicate masked partition')
        output = np.full(len(blocks), np.nan)
        records = []
        folds[key] = dict(profile=output, records=records)
        def capture(lo, hi, costs, medians, unit, denominator, penalty):
            prediction, record = _average(costs, medians, unit, denominator, penalty)
            assert np.isnan(output[lo:hi]).all()
            output[lo:hi] = prediction
            records.append(record)
        return capture
    captured = dict(partition=partition)
    baseline, old_decision, _ = gene_segments._select(values, gene_bins, blocks, capture=captured, solver=solver)
    if 'P' not in captured:
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
    weights, record = _convex_l1(P, target, solver=solver)
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



