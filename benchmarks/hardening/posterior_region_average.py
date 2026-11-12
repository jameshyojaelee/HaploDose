"""Exact averaging of descriptive RNA region summaries over HMM boundaries."""
import argparse
from functools import lru_cache
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from haplodose import profile
from haplodose.io import write_json


@lru_cache(None)
def load(name):
    spec = importlib.util.spec_from_file_location('posterior_regions_'+name,
        Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_probabilities(emission, prior, change):
    """Probability each interval is one complete constant-state run."""
    e, prior = np.asarray(emission, float), np.asarray(prior, float)
    if (e.ndim != 2 or not len(e) or prior.shape != (e.shape[1],)
            or not np.isfinite(e).all() or not np.isfinite(prior).all()
            or np.any(prior <= 0) or not 0 < change < 1):
        raise ValueError('Finite emissions, positive prior and interior change probability required')
    prior = prior/prior.sum()
    n, states = e.shape
    T = np.log(load('shared_membership_hmm').transition(prior, change))
    forward, backward = np.empty_like(e), np.empty_like(e)
    forward[0] = np.log(prior)+e[0]
    for j in range(1, n):
        forward[j] = e[j]+logsumexp(forward[j-1, :, None]+T, axis=0)
    backward[-1] = 0.
    for j in range(n-2, -1, -1):
        backward[j] = logsumexp(T+e[j+1, None, :]+backward[j+1, None, :], axis=1)
    evidence = float(logsumexp(forward[-1]))
    off = T.copy()
    np.fill_diagonal(off, -np.inf)
    left, right = np.empty_like(e), np.empty_like(e)
    left[0], right[-1] = np.log(prior), 0.
    for j in range(1, n):
        left[j] = logsumexp(forward[j-1, :, None]+off, axis=0)
    for j in range(n-1):
        right[j] = logsumexp(off+e[j+1, None, :]+backward[j+1, None, :], axis=1)
    cumulative = np.vstack([np.zeros(states), np.cumsum(e, axis=0)])
    probabilities = np.zeros((n, n))
    for a in range(n):
        for b in range(a, n):
            logp = left[a]+cumulative[b+1]-cumulative[a]+(b-a)*np.diag(T)+right[b]-evidence
            probabilities[a, b] = float(np.exp(logsumexp(logp)))
    mass = np.array([math.fsum(map(float, probabilities[:j+1, j:].ravel())) for j in range(n)])
    # Bound accumulated log-operation error, its exp propagation and mass sum.
    operations = 128*n*n*states+256
    gamma = operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)
    scale = max(1., abs(evidence), float(np.max(abs(cumulative))), float(np.max(abs(forward))))
    budget = math.expm1(8*gamma*scale)+gamma+n*n*np.nextafter(0., 1.)
    error = float(np.max(abs(mass-1)))
    if not math.isfinite(budget) or error > budget:
        raise ArithmeticError('Unresolved posterior run-mass conservation')
    return probabilities, dict(bins=n, states=states, maximum_mass_error=error,
        mass_rounding_budget=budget, log_evidence=evidence)


class MassTree:
    """Integer dyadic masses: exact weighted-median ties without float epsilons."""
    def __init__(self, n):
        self.tree = [0]*(n+1)
        self.total = 0

    def add(self, position, value):
        self.total += value
        j = position+1
        while j < len(self.tree):
            self.tree[j] += value
            j += j & -j

    def locate(self, target):
        index, cumulative = 0, 0
        step = 1 << (len(self.tree).bit_length()-1)
        while step:
            nxt = index+step
            if nxt < len(self.tree) and cumulative+self.tree[nxt] < target:
                cumulative += self.tree[nxt]
                index = nxt
            step >>= 1
        return index, cumulative

    def median(self, ordered):
        index, _ = self.locate((self.total+1)//2)
        # The first positive mass at/after index is recovered exactly.
        after = 0
        j = index+1
        while j:
            after += self.tree[j]
            j -= j & -j
        if 2*after == self.total:
            following, _ = self.locate(after+1)
            return float((ordered[index]+ordered[following])/2)
        return float(ordered[index])


def interval_summaries(values, variance, bins, length):
    values, variance, bins = np.asarray(values, float), np.asarray(variance, float), np.asarray(bins, int)
    if (values.ndim != 1 or variance.shape != values.shape or bins.shape != values.shape
            or not np.isfinite(values).all() or not np.isfinite(variance).all() or np.any(variance <= 0)
            or not np.array_equal(np.unique(bins), np.arange(length))):
        raise ValueError('Finite gene effects, positive weights and complete local bins required')
    order = np.argsort(values, kind='stable')
    y, v, b = values[order], variance[order], bins[order]
    weights = 1/np.sqrt(v)
    dyadic = [float(w).as_integer_ratio() for w in weights]
    exponent = max(d.bit_length()-1 for _, d in dyadic)
    integers = [p << (exponent-(q.bit_length()-1)) for p, q in dyadic]
    gene_rows = [np.flatnonzero(b == j) for j in range(length)]
    summaries = np.zeros((3, length, length))
    for a in range(length):
        weighted, uniform = MassTree(len(y)), MassTree(len(y))
        selected = []
        for end in range(a, length):
            for j in gene_rows[end]:
                weighted.add(int(j), integers[j])
                uniform.add(int(j), 1)
                selected.append(int(j))
            summaries[0, a, end] = uniform.median(y)
            summaries[1, a, end] = weighted.median(y)
            summaries[2, a, end] = math.fsum(float(y[j]/v[j]) for j in selected)/math.fsum(float(1/v[j]) for j in selected)
    return summaries


def average(values, variance, gene_bins, reset, fitted, params):
    y, total = np.asarray(fitted['bin_values']), np.asarray(fitted['total_variance'])
    means = fitted['beta']*np.asarray(profile.LOG_DOSE, float)
    emission = -.5*(np.log(2*np.pi*total)[:, None]+(y[:, None]-means)**2/total[:, None])
    starts = np.flatnonzero(reset)
    ends = np.r_[starts[1:], len(y)]
    result = np.empty((3, len(y)))
    records = []
    for a, b in zip(starts, ends):
        probability, record = run_probabilities(emission[a:b], profile.STATE_PRIOR, params.p_change)
        keep = (gene_bins >= a)&(gene_bins < b)
        summaries = interval_summaries(values[keep], variance[keep], gene_bins[keep]-a, b-a)
        for j in range(b-a):
            for k in range(3):
                result[k, a+j] = math.fsum(map(float,
                    (probability[:j+1, j:]*summaries[k, :j+1, j:]).ravel()))
        record.update(first_bin=int(a), last_bin=int(b-1))
        records.append(record)
    return result, records


def checks():
    gaussian = load('reference_gaussian_regions')
    prior, change = np.array([.2, .5, .3]), .1
    emission = np.log(np.array([[.8, .2, .5], [.1, .6, .3], [.2, .5, .7]]))
    probability, record = run_probabilities(emission, prior, change)
    T = load('shared_membership_hmm').transition(prior, change)
    paths = list(itertools.product(range(3), repeat=3))
    weights = []
    for path in paths:
        value = prior[path[0]]
        for j, s in enumerate(path):
            value *= math.exp(emission[j, s])
            if j:value *= T[path[j-1], s]
        weights.append(value)
    weights = np.asarray(weights)/math.fsum(weights)
    independent = np.zeros((3, 3))
    y, v, bins = np.array([0., 2., 1., 4., -1., 5.]), np.array([1., 4., 1., 9., 4., 1.]), np.repeat(np.arange(3), 2)
    summaries = interval_summaries(y, v, bins, 3)
    exact = load('reference_weighted_regions')
    reference = np.zeros((3, 3))
    for path, weight in zip(paths, weights):
        for a in range(3):
            if a and path[a] == path[a-1]:continue
            end = a
            while end+1 < 3 and path[end+1] == path[a]:end += 1
            independent[a, end] += weight
            keep = (bins >= a)&(bins <= end)
            values = [float(np.median(y[keep])), exact.weighted_median(y[keep], 1/np.sqrt(v[keep])),
                math.fsum(map(float, y[keep]/v[keep]))/math.fsum(map(float, 1/v[keep]))]
            np.testing.assert_array_equal(summaries[:, a, end], values)
            for j in range(a, end+1):reference[:, j] += weight*np.asarray(values)
    error = float(np.max(abs(probability-independent)))
    assert error <= record['mass_rounding_budget']
    direct = np.array([[math.fsum(map(float, (probability[:j+1,j:]*summaries[k,:j+1,j:]).ravel()))
        for j in range(3)] for k in range(3)])
    mean_error = float(np.max(abs(direct-reference)))
    assert mean_error <= record['mass_rounding_budget']*max(1., float(np.max(abs(y))))
    # A one-bin block has exactly one run, so its summaries are unchanged.
    one, one_record = run_probabilities(emission[:1], prior, change)
    assert abs(one[0,0]-1) <= one_record['mass_rounding_budget']
    alternate = interval_summaries(y[::-1], v[::-1], bins[::-1], 3)
    np.testing.assert_array_equal(summaries, alternate)
    return dict(exhaustive_paths=len(paths), maximum_run_probability_error=error,
        maximum_averaged_summary_error=mean_error, exact_interval_medians=True,
        cell_state_probabilities_calibrated=False, gene_permutation_exact=True,
        one_bin_limit=True, mass=record, original_gaussian_checks=gaussian.checks())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--checks-only', action='store_true')
    parser.add_argument('--reference-refit', type=int, choices=range(9), default=0)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    if args.checks_only:write_json(checks(), args.out/'report.json')
    else:
        weighted = load('reference_weighted_regions')
        comparison = args.task.parent/'normal_noise_decomposition_20261109'
        if args.reference_refit:
            # Carry the full39-method family through each existing reference
            # resample; never drop earlier comparators for missing refits.
            centered = args.task/'resampled_centered'/f'refit_{args.reference_refit:03d}_v1'
            centered.mkdir(parents=True, exist_ok=False)
            weighted.panel(argparse.Namespace(out=centered, task=args.task,
                reference_refit=args.reference_refit,
                comparison_root=args.task.parent/'reference_gaussian_regions_20261109',
                weighting='reference', normal_centering=True,
                reference_gaussian_regions=True, gaussian_centre='bin-median'))
            decomposed = args.task/'resampled_decomposed'/f'refit_{args.reference_refit:03d}_v1'
            decomposed.mkdir(parents=True, exist_ok=False)
            weighted.panel(argparse.Namespace(out=decomposed, task=args.task,
                reference_refit=args.reference_refit, comparison_root=args.task/'resampled_centered',
                weighting='reference', normal_centering=True,
                reference_gaussian_regions=True, gaussian_variance='noise-subtracted'))
            comparison = args.task/'resampled_decomposed'
        weighted.panel(argparse.Namespace(out=args.out, task=args.task,
            reference_refit=args.reference_refit, comparison_root=comparison,
            weighting='reference', normal_centering=True, reference_gaussian_regions=True,
            posterior_boundaries=True))
        if args.reference_refit == 0 and args.out.name != 'bulk_v1':
            original = args.task/'bulk_v1'
            with np.load(original/'RNA_predictions.npz') as old, np.load(args.out/'RNA_predictions.npz') as new:
                assert old.files == new.files
                for key in old.files:np.testing.assert_array_equal(old[key], new[key])
            pd.testing.assert_frame_equal(pd.read_csv(original/'patient_scores.tsv', sep='\t', float_precision='round_trip'),
                pd.read_csv(args.out/'patient_scores.tsv', sep='\t', float_precision='round_trip'), check_exact=True)
            write_json(dict(primary_predictions_and_all39_scores_exact=True), args.out/'refit_control.json')


if __name__ == '__main__':main()
