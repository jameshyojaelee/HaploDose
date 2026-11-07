"""Prespecified oracle-model simultaneous-interval experiment, compute nodes only."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import special, stats

from haplodose.profile import ProfileParams, fit_pooled, segment, segments_table
from haplodose.provenance import start_run, finish_run


def file_record(path):
    p = Path(path).resolve()
    return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def exact(successes, total, confidence=.95):
    ci = stats.binomtest(int(successes), int(total)).proportion_ci(confidence_level=confidence, method='exact')
    return dict(successes=int(successes), n=int(total), fraction=successes / total,
                cp_low=float(ci.low), cp_high=float(ci.high), confidence=confidence)


def coverage_family(rows, target, alpha):
    """Original six-member exact undercoverage test, retaining extreme log tails."""
    if len(rows) != 6:
        raise RuntimeError('Incomplete prespecified six-endpoint family')
    logp = np.array([special.logsumexp(stats.binom.logpmf(
        np.arange(row['coverage']['successes'] + 1), row['coverage']['n'], target)) for row in rows])
    order = np.argsort(logp, kind='stable')
    adjusted = np.empty(6)
    adjusted[order] = np.minimum(0., np.maximum.accumulate(logp[order] + np.log(np.arange(6, 0, -1))))
    for row, lp, la in zip(rows, logp, adjusted):
        row.update(log_p_undercoverage=float(lp), log_p_holm=float(la),
                   p_undercoverage=float(np.exp(lp)), p_holm=float(np.exp(la)),
                   falsified_undercoverage=bool(la < np.log(alpha)),
                   simultaneous_cp=exact(row['coverage']['successes'], row['coverage']['n'], 1 - alpha / 6))


def summarize(frame, expected_n):
    if len(frame) != expected_n or frame[['covered', 'excludes_zero', 'width']].isna().any().any():
        raise RuntimeError('Missing required interval outcomes')
    return dict(coverage=exact(frame.covered.sum(), len(frame)),
                zero_exclusion=exact(frame.excludes_zero.sum(), len(frame)),
                infinite_width=exact(np.isinf(frame.width).sum(), len(frame)),
                median_width=float(np.median(frame.width)),
                finite_width_n=int(np.isfinite(frame.width).sum()))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--generator-source', type=Path, required=True)
    parser.add_argument('--design', type=Path, required=True)
    args = parser.parse_args(argv)
    design = json.loads(args.design.read_text())
    if file_record(args.generator_source)['sha256'] != design['generator_sha256']:
        raise ValueError('Original independent generator source hash does not match frozen design')
    module_spec = importlib.util.spec_from_file_location('original_uncertainty', args.generator_source)
    original = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(original)
    out = original.fresh_output(args.task_root, args.out)
    provenance = start_run(out, {'script': __file__, 'generator_source': str(args.generator_source),
                                 'design': str(args.design)}, design,
                           {'measurement': 'oracle Gaussian log-intensity plus Poisson counts',
                            'replication': 'simulation draws, no biological replication'}, mode='uncertainty_simultaneous')
    try:
        import jax
        from haplodose.device import select
        backend = select('cpu')
        devices = [dict(platform=d.platform, kind=d.device_kind) for d in jax.devices()]
        if backend != 'cpu' or any(d['platform'] != 'cpu' for d in devices):
            raise RuntimeError('Uncertainty experiment requires the recorded CPU inference backend')
        (out / 'backend.json').write_text(json.dumps(dict(backend=backend, devices=devices), indent=2))
        experiment(original, design, out)
        (out / 'source_manifest.json').write_text(json.dumps(dict(
            executed_script=file_record(__file__), design=file_record(args.design),
            generator=file_record(args.generator_source), command=sys.argv), indent=2))
        finish_run(out, provenance, backend=backend)
    except Exception as error:
        finish_run(out, provenance, error=error)
        raise


def experiment(original, design, out):
    from haplodose.uncertainty import systematic_rna_intervals
    J, fixed = design['bins'], design['fixed_bins']
    per_chrom = J // 4
    reset = np.zeros(J, bool)
    reset[::per_chrom] = True
    params = ProfileParams()
    fixed_a = (per_chrom - fixed) // 2
    heterogeneous = np.tile(np.r_[np.full(per_chrom // 2, .25 ** 2),
                                 np.full(per_chrom - per_chrom // 2, .02 ** 2)], 4)
    settings = [('matched_homogeneous', np.zeros(J), 0.),
                ('matched_heterogeneous', heterogeneous, 0.),
                ('correlated_heterogeneous', heterogeneous, .8)]
    if [s[0] for s in settings] != design['conditions']:
        raise ValueError('Prespecified conditions changed')
    tau_true, records, saved = .08, [], {}
    z = float(stats.norm.ppf((1 + design['coverage_target']) / 2))
    effects = [float(m * params.beta * np.log(1.5)) for m in design['power_effect_multipliers']]
    for phase in ('coverage', 'power'):
        n = design['replicates_per_condition'] if phase == 'coverage' else design['power_replicates_per_condition_effect']
        for ci, (condition, reference_var, rho) in enumerate(settings):
            for ei, effect in enumerate([0.] if phase == 'coverage' else effects):
                count_store, latent_store = [], []
                for replicate in range(n):
                    seed = [design['seed'], ci, replicate] if phase == 'coverage' else [design['seed'], 10, ci, replicate]
                    rng = np.random.default_rng(seed)
                    y, latent, expected, conditional_mean, covariance = original.generate(rng, J, tau_true, reference_var, rho)
                    mu = np.zeros(J)
                    mu[fixed_a:fixed_a + fixed] = effect
                    if effect:
                        conditional_mean = conditional_mean * np.exp(mu)
                        y = rng.poisson(conditional_mean).astype(float)
                    center, beta, tau = fit_pooled(y, expected, reset, params, v=reference_var)
                    post, path, dose, _ = segment(y[None], expected[None], reset, [center], beta,
                                                  np.sqrt(tau ** 2 + reference_var)[None], params)
                    reported = segments_table(post[0], path[0], dose[0], y, expected, reset, beta, tau, params,
                                              reference_var=reference_var)
                    if not reported:
                        raise RuntimeError('Missing actual HMM runs')
                    observed = np.log((y + .5) / (expected + .5))
                    winner = max(reported, key=lambda r: abs(observed[r['first_bin']:r['last_bin'] + 1].mean()))
                    intervals = [(fixed_a, fixed_a + fixed - 1), (winner['first_bin'], winner['last_bin'])]
                    # API is explicit: oracle covariance and independent expected
                    # median intensity, not an estimated-reference plug-in.
                    bounds = simultaneous_bounds(systematic_rna_intervals, y, expected, covariance, reset, intervals,
                                                 1 - design['coverage_target'])
                    for selection, (a, b), (lower, upper) in zip(design['intervals'], intervals, bounds):
                        target = float(mu[a:b + 1].mean())
                        old = original.interval_record(y, latent, conditional_mean, covariance, dose[0], a, b,
                                                       tau, reference_var, beta, params, z)
                        estimate, se = old['observed_mean_rna_log_ratio'], old['reported_conditional_se']
                        for method, lo, hi in [('legacy_reported_se', estimate - z * se, estimate + z * se),
                                               ('simultaneous', lower, upper)]:
                            if np.isnan([lo, hi]).any() or lo > hi or not np.isfinite([target, estimate, se]).all():
                                raise RuntimeError('Missing or invalid interval output')
                            records.append(dict(phase=phase, condition=condition, effect=effect, replicate=replicate,
                                interval=selection, method=method, first_bin=a, last_bin=b, true_mu_I=target,
                                lower=float(lo), upper=float(hi), width=float(hi - lo),
                                covered=bool(lo <= target <= hi), excludes_zero=bool(lo > 0 or hi < 0),
                                estimated_tau=float(tau), estimated_beta=float(beta)))
                    count_store.append(y)
                    latent_store.append(latent)
                key = f'{phase}_{condition}_effect{ei}'
                saved[key + '_counts'] = np.asarray(count_store)
                saved[key + '_latent'] = np.asarray(latent_store)
                saved[key + '_covariance'] = covariance
                saved[key + '_systematic_mu'] = mu
                print(json.dumps(dict(phase=phase, condition=condition, effect=effect,
                                      simulations_complete=n)), flush=True)
    frame = pd.DataFrame(records)
    frame.to_csv(out / 'intervals.tsv', sep='\t', index=False)
    np.savez_compressed(out / 'simulations.npz', **saved)
    families = {}
    for method in ('legacy_reported_se', 'simultaneous'):
        rows = []
        for condition in design['conditions']:
            for interval in design['intervals']:
                subset = frame[(frame.phase == 'coverage') & frame.method.eq(method) & frame.condition.eq(condition)
                               & frame.interval.eq(interval)]
                rows.append(dict(condition=condition, interval=interval, **summarize(subset, design['replicates_per_condition'])))
        coverage_family(rows, design['coverage_target'], design['family_alpha'])
        families[method] = rows
    power = []
    for (condition, effect, interval, method), subset in frame[frame.phase == 'power'].groupby(
            ['condition', 'effect', 'interval', 'method'], sort=False):
        power.append(dict(condition=condition, effect=float(effect), interval=interval, method=method,
                          **summarize(subset, design['power_replicates_per_condition_effect'])))
    (out / 'summary.json').write_text(json.dumps(dict(design=design, profile_params=asdict(params),
        primary_family='simultaneous', separate_comparator_family='legacy_reported_se',
        coverage_families=families, secondary_operating_curves=power, effect_grid=effects,
        limitation='Oracle supplied covariance and independent median intensity only; no reference-fit or biological validation. '
                   'Nonrejection of undercoverage does not certify coverage; power curves are descriptive with no universal threshold.'), indent=2))


def simultaneous_bounds(helper, y, expected, covariance, reset, intervals, alpha):
    """Use generating covariance, not a fitted tau/reference substitution."""
    result = helper(y, expected, covariance, reset, intervals,
                    counts_kind='raw', count_provenance='Untransformed integer Poisson draws from frozen generator; no pooling or normalization',
                    expectation_provenance='Frozen original generator deterministic median intensity E=2000 independent of observations',
                    covariance_provenance='Exact covariance returned by frozen independent Gaussian log-intensity generator',
                    bin_universe_provenance='Prespecified 120 bins in four fixed 30-bin chromosomes; design frozen before observations',
                    confidence=1 - alpha)
    rows = result['intervals']
    if len(rows) != len(intervals):
        raise RuntimeError('Missing simultaneous helper interval outputs')
    for row, (a, b) in zip(rows, intervals):
        if row['first_bin'] != a or row['last_bin'] != b:
            raise RuntimeError('Helper interval alignment differs from frozen selection')
    return [(row['lower'], row['upper']) for row in rows]


if __name__ == '__main__':
    main()
