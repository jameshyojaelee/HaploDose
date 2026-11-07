"""Known-variance fixed/selected-interval falsifiers; run on a compute node.

Independent simulations measure conditional Monte Carlo error, not biological
replication. Replicate count is explicit; exact binomial intervals report the
achieved precision instead of treating a small pilot as a calibration gate.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np
from scipy import special, stats

from haplodose.io import new_output_dir, write_json, write_table
from haplodose.profile import ProfileParams, fit_pooled, segment, segments_table
from haplodose.provenance import finish_run, start_run


def fresh_output(task_root, output):
    root, output = task_root.resolve(strict=True), output.resolve()
    relative = output.relative_to(root)
    if not relative.parts or relative.parts[0] == 'baseline':
        raise ValueError('Probe output must be a new descendant of task_root outside baseline')
    if output.exists():
        raise FileExistsError(output)
    return new_output_dir(output)


def binomial_summary(successes, total, confidence):
    interval = stats.binomtest(successes, total).proportion_ci(confidence_level=confidence, method='exact')
    return dict(successes=int(successes), n_simulations=int(total), fraction=successes / total,
                exact_ci_low=float(interval.low), exact_ci_high=float(interval.high), confidence=confidence,
                interval_scope='nominal_descriptive_not_simultaneous')


def generate(rng, J, tau, reference_var, rho):
    """Independent generator: Gaussian log-intensity, optionally correlated by chromosome.

    E is a median intensity, matching the documented PLN parameterization. It is
    not silently interpreted as the marginal arithmetic mean of the PLN.
    """
    covariance = np.zeros((J, J))
    per_chrom = J // 4
    sd = np.sqrt(tau ** 2 + reference_var)
    for a in range(0, J, per_chrom):
        b = a + per_chrom
        distance = np.abs(np.arange(per_chrom)[:, None] - np.arange(per_chrom)[None])
        covariance[a:b, a:b] = sd[a:b, None] * sd[None, a:b] * rho ** distance
    latent = rng.multivariate_normal(np.zeros(J), covariance)
    expected = np.full(J, 2000.)
    mean = expected * np.exp(latent)
    return rng.poisson(mean).astype(float), latent, expected, mean, covariance


def interval_record(y, latent, conditional_mean, covariance, dose, a, b, estimated_tau,
                    reference_var, beta, params, z):
    sl = slice(a, b + 1)
    n = b - a + 1
    post = np.zeros((n, 4))
    post[:, 1] = 1
    rows = segments_table(post, np.ones(n, int), dose[sl], y[sl], np.full(n, 2000.),
                          np.r_[True, np.zeros(n - 1, bool)], beta, estimated_tau, params,
                          reference_var=reference_var[sl])
    observed = np.log((y[sl] + .5) / (2000. + .5)).mean()
    se = rows[0]['se_rna_log_ratio_approx']
    # For latent Gaussian offsets this oracle variance is exact; the additional
    # conditional Poisson log-count variance remains a delta approximation.
    latent_se = np.sqrt(covariance[sl, sl].sum()) / n
    count_se = np.sqrt(covariance[sl, sl].sum() + np.sum(1 / conditional_mean[sl])) / n
    return dict(first_bin=a, last_bin=b, n_bins=n, observed_mean_rna_log_ratio=float(observed),
                posterior_mean_rna_effect=float(dose[sl].mean()), reported_conditional_se=se,
                oracle_latent_se=float(latent_se), oracle_log_count_se_approx=float(count_se),
                latent_mean=float(latent[sl].mean()),
                reported_mean_covers_zero=bool(abs(observed) <= z * se),
                reported_posterior_covers_zero=bool(abs(dose[sl].mean()) <= z * se),
                oracle_latent_covers_zero=bool(abs(latent[sl].mean()) <= z * latent_se),
                oracle_log_count_covers_zero=bool(abs(observed) <= z * count_se))


def run(args, out):
    import pandas as pd

    J, per_chrom = args.bins, args.bins // 4
    reset = np.zeros(J, bool)
    reset[::per_chrom] = True
    params = ProfileParams()
    tau_true = .08
    settings = [('matched_homogeneous', np.zeros(J), 0.),
                ('matched_heterogeneous', np.tile(np.r_[np.full(per_chrom // 2, .25 ** 2),
                    np.full(per_chrom - per_chrom // 2, .02 ** 2)], 4), 0.),
                ('correlated_heterogeneous', np.tile(np.r_[np.full(per_chrom // 2, .25 ** 2),
                    np.full(per_chrom - per_chrom // 2, .02 ** 2)], 4), .8)]
    records, arrays = [], {}
    z = float(stats.norm.ppf((1 + args.confidence) / 2))
    for condition_index, (condition, reference_var, rho) in enumerate(settings):
        count_store, latent_store = [], []
        for replicate in range(args.replicates):
            rng = np.random.default_rng([args.seed, condition_index, replicate])
            y, latent, expected, conditional_mean, covariance = generate(rng, J, tau_true, reference_var, rho)
            center, beta, tau = fit_pooled(y, expected, reset, params, v=reference_var)
            post, path, dose, _ = segment(y[None], expected[None], reset, [center], beta,
                                          np.sqrt(tau ** 2 + reference_var)[None], params)
            reported = segments_table(post[0], path[0], dose[0], y, expected, reset, beta, tau, params,
                                      reference_var=reference_var)
            fixed_a = (per_chrom - args.fixed_bins) // 2
            fixed = interval_record(y, latent, conditional_mean, covariance, dose[0], fixed_a,
                                    fixed_a + args.fixed_bins - 1, tau, reference_var, beta, params, z)
            # Explicit winner selection over the actual reported HMM runs. Neutral
            # runs remain eligible; this is not disguised as a discovery threshold.
            observed_r = np.log((y + .5) / (expected + .5))
            selected_run = max(reported, key=lambda r: abs(observed_r[r['first_bin']:r['last_bin'] + 1].mean()))
            selected = interval_record(y, latent, conditional_mean, covariance, dose[0],
                                       selected_run['first_bin'], selected_run['last_bin'], tau,
                                       reference_var, beta, params, z)
            base = dict(condition=condition, replicate=replicate, true_tau=tau_true, estimated_tau=tau,
                        rho=rho, center=center, beta=beta,
                        any_resolved_event=any(row['resolved'] for row in reported),
                        variance_mode='conservative_roughness_plus_reference')
            records.extend([{**base, 'interval': 'prespecified', **fixed},
                            {**base, 'interval': 'largest_absolute_mean_among_HMM_runs', **selected}])
            count_store.append(y)
            latent_store.append(latent)
        arrays[condition + '_counts'] = np.vstack(count_store)
        arrays[condition + '_latent'] = np.vstack(latent_store)
        arrays[condition + '_reference_var'] = reference_var
    write_table(pd.DataFrame(records), out / 'intervals.tsv')
    np.savez_compressed(out / 'simulations.npz', **arrays)
    summary = []
    for condition, _, _ in settings:
        for interval in ('prespecified', 'largest_absolute_mean_among_HMM_runs'):
            rows = [r for r in records if r['condition'] == condition and r['interval'] == interval]
            row = dict(condition=condition, interval=interval,
                       mean_estimated_tau=float(np.mean([r['estimated_tau'] for r in rows])),
                       mean_reported_se=float(np.mean([r['reported_conditional_se'] for r in rows])),
                       mean_oracle_se_approx=float(np.mean([r['oracle_log_count_se_approx'] for r in rows])))
            for key in ('reported_mean_covers_zero', 'reported_posterior_covers_zero',
                        'oracle_latent_covers_zero', 'oracle_log_count_covers_zero', 'any_resolved_event'):
                row[key] = binomial_summary(sum(r[key] for r in rows), len(rows), args.confidence)
            summary.append(row)
    # Predeclared primary family: this endpoint only, across all three generators
    # and both fixed/selected intervals. The remaining endpoints stay descriptive.
    family_alpha = .05
    family_size = 6
    if len(summary) != family_size:
        raise RuntimeError('The predeclared six-member undercoverage family is incomplete')
    # Binomial tails can be smaller than float64's smallest positive value.
    # Retain log probabilities so underflow is never interpreted as p=0.
    log_p = np.array([special.logsumexp(stats.binom.logpmf(
        np.arange(row['reported_mean_covers_zero']['successes'] + 1), args.replicates, args.confidence))
        for row in summary])
    order = np.argsort(log_p, kind='stable')
    log_adjusted = np.empty(family_size)
    log_adjusted[order] = np.minimum(0., np.maximum.accumulate(log_p[order] + np.log(np.arange(family_size, 0, -1))))
    for row, log_probability, log_holm in zip(summary, log_p, log_adjusted):
        primary = row['reported_mean_covers_zero']
        interval = stats.binomtest(primary['successes'], args.replicates).proportion_ci(
            confidence_level=1 - family_alpha / family_size, method='exact')
        primary.update(p_undercoverage_exact=float(np.exp(log_probability)), p_undercoverage_holm=float(np.exp(log_holm)),
                       log_p_undercoverage_exact=float(log_probability), log_p_undercoverage_holm=float(log_holm),
                       p_value_underflow=bool(log_probability < np.log(np.nextafter(0., 1.))),
                       falsified_undercoverage=bool(log_holm < np.log(family_alpha)),
                       simultaneous_cp_ci_low=float(interval.low), simultaneous_cp_ci_high=float(interval.high),
                       simultaneous_family_confidence=1 - family_alpha,
                       simultaneous_interval_method='two-sided Clopper-Pearson; Bonferroni over six primary coverages')
    write_json(dict(summary=summary, profile_params=asdict(params), true_tau=tau_true,
                    replicates=args.replicates, seed=args.seed, confidence=args.confidence,
                    fixed_interval_bins=args.fixed_bins, bins=J,
                    oracle_definition='Exact Gaussian latent-offset variance; Poisson log-count variance uses a delta approximation',
                    selected_interval_definition='Largest absolute observed mean among actual HMM runs, including neutral runs',
                    replication='Independent simulated specimens conditional on generator; no biological replication',
                    primary_undercoverage_family=dict(endpoint='reported_mean_covers_zero', size=family_size,
                        conditions=[condition for condition, _, _ in settings],
                        intervals=['prespecified', 'largest_absolute_mean_among_HMM_runs'],
                        nominal_coverage_target=args.confidence, family_alpha=family_alpha,
                        test='one-sided exact binomial undercoverage test against nominal coverage target',
                        multiplicity='Holm over all six primary tests',
                        intervals_method='two-sided Clopper-Pearson, Bonferroni simultaneous family confidence 0.95',
                        decision='falsified_undercoverage when Holm p < 0.05; no automatic scientific promotion',
                        other_endpoints='nominal descriptive intervals; excluded from the primary test family'),
                    calibration_status='not_established',
                    limitations=['Numerical/variance diagnostic, not a benchmark of biological CN accuracy',
                                 'No independently estimated reference panel or programme parameters',
                                 'Conservative tau/v mode retained; no automatic accuracy claim from an algebraic change',
                                 'Primary six-member family can falsify synthetic undercoverage; nonrejection does not certify error control',
                                 'Correlated offsets are an explicit misspecification scenario, not a claim about a particular tissue']),
               out / 'summary.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--replicates', type=int, required=True)
    parser.add_argument('--seed', type=int, default=20261106)
    parser.add_argument('--bins', type=int, default=120)
    parser.add_argument('--fixed-bins', type=int, default=6)
    parser.add_argument('--confidence', type=float, default=.95)
    args = parser.parse_args()
    if (args.replicates < 1 or args.bins < 12 or args.bins % 4 or not 1 <= args.fixed_bins <= args.bins // 4
            or not 0 < args.confidence < 1):
        parser.error('require replicates>=1, bins>=12 divisible by four, 1<=fixed-bins<=bins/4, confidence in (0,1)')
    out = fresh_output(args.task_root, args.out)
    record = start_run(out, {'probe_source': __file__},
                       {k: v for k, v in vars(args).items() if k not in ('out', 'task_root')},
                       {'measurement': 'independently generated Poisson-lognormal counts',
                        'replication': 'conditional Monte Carlo'}, mode='uncertainty_diagnostic')
    try:
        run(args, out)
        finish_run(out, record, backend='cpu')
        print(json.dumps({'summary': str(out / 'summary.json'), 'replicates_per_condition': args.replicates}))
    except Exception as error:
        finish_run(out, record, error=error)
        raise


if __name__ == '__main__':
    main()
