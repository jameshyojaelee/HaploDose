"""Falsification diagnostic: selected-parent filter versus complete workflow replay.

Run on a compute node, e.g. python benchmarks/hardening/clones.py --out NEW_DIRECTORY.
Three conditions retain a truncal RNA dose: matched NB null, continuous correlated state (one DNA clone),
and two DNA clones. Independent simulations, not patients, are the replication units. This is a small
expression-only diagnostic; it does not certify a biological false-split rate or calibrate the allele channel.
Complete replay refits the entire clone workflow, conditional on a known reference and fixed hyperparameters.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import hashlib
from importlib import metadata
import json
from pathlib import Path
import platform
import sys

import numpy as np
from scipy import stats
from sklearn.metrics import adjusted_rand_score

from haplodose.clones import CloneParams, Rows, SpecimenModel, depth_moments, fit_specimen

CONDITIONS = ('matched_nb', 'correlated_state', 'two_clones')
METHODS = ('raw', 'selected_filter', 'full_replay')
GENERATOR_VERSION = 'same_nuisance_v2'
NUISANCE_MODEL = dict(region='bins J//2:3*J//4', distribution='independent per-cell Normal(0, 0.25**2)',
                      normalization='rates renormalized across all bins after state and CNA effects',
                      applies_to=['correlated_state', 'two_clones'])


def array_hash(array):
    array = np.ascontiguousarray(array)
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode())
    digest.update(json.dumps(array.shape).encode())
    digest.update(array.tobytes())
    return digest.hexdigest()


def generate(rng, condition, n, J, *, return_latent=False):
    base = np.full(J, 1. / J)
    phi = np.full(J, .15)
    depth = np.exp(rng.normal(np.log(800.), .3, n))
    dose = np.zeros((n, J))
    dose[:, :J // 4] = .6 * np.log(1.5)
    truth = np.zeros(n, int)
    state_effect = np.zeros(n)
    if condition in ('correlated_state', 'two_clones'):
        state_effect = rng.normal(0., .25, n)
        dose[:, J // 2:3 * J // 4] += state_effect[:, None]
    if condition == 'two_clones':
        truth[n // 2:] = 1
        dose[n // 2:, J // 2:3 * J // 4] += .6 * np.log(1.5)
    elif condition not in ('matched_nb', 'correlated_state'):
        raise ValueError(condition)
    rate = base[None] * np.exp(dose)
    rate /= rate.sum(axis=1, keepdims=True)
    mean = depth[:, None] * rate
    Y = rng.negative_binomial(1 / phi, 1 / (1 + phi * mean)).astype(np.float32)
    output = (Y, base, phi, truth)
    return output + (dict(mean=mean, log_rate_effect=dose, depth=depth, state_effect=state_effect),) if return_latent else output


def fit(Y, base, phi, reset, params):
    # Counts are the actual modeled bin counts; no synthetic gene-veto claim is made.
    return fit_specimen(Rows(Y, np.arange(len(Y))), Y.sum(axis=1, dtype=np.float64), base[:, None],
                        reset, phi, params)


def statistic(result, min_cells):
    """Number of non-reference MAP groups reaching min_cells, not number of validated DNA clones.

    Soft responsibilities below MAP size do not count. More than one group in the two-clone alternative
    measures multigroup detection only: over-splitting and incorrect cell membership can still satisfy it.
    """
    labels = result['R'].argmax(axis=1)
    sizes = np.bincount(labels, minlength=len(result['doses']))
    return int((sizes[1:] >= min_cells).sum())


def bounds(successes, n, alpha=.05):
    """One-sided exact binomial lower and upper bounds, each at confidence 1-alpha."""
    if not 0 <= successes <= n or n < 1 or not 0 < alpha < 1:
        raise ValueError('need 0 <= successes <= n, n >= 1 and 0 < alpha < 1')
    lower = 0. if successes == 0 else float(stats.beta.ppf(alpha, successes, n - successes + 1))
    upper = 1. if successes == n else float(stats.beta.ppf(1 - alpha, successes + 1, n - successes))
    return lower, upper


def self_check():
    """Small mathematical/endpoint invariants; no model fitting."""
    result = dict(R=np.eye(3)[[0, 0, 1, 1, 2]], doses=np.zeros((3, 2)))
    assert statistic(result, 2) == 1       # reference group excluded; one-cell group excluded
    assert statistic(result, 1) == 2
    assert bounds(0, 59)[1] <= .05
    assert bounds(0, 59, .05 / 6)[1] > .05  # single-endpoint precision is insufficient for six null endpoints
    assert bounds(0, 94, .05 / 6)[1] <= .05
    assert bounds(59, 59)[0] > .9
    assert bounds(0, 10)[0] == 0. and bounds(10, 10)[1] == 1.


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--replicates', type=int, required=True,
                        help='Declare simulation count before execution; small counts are screening only.')
    parser.add_argument('--condition', choices=('all',) + CONDITIONS, default='all',
                        help='One condition permits independent ordinary jobs; no array is needed.')
    parser.add_argument('--draws', type=int, default=19)
    parser.add_argument('--cells', type=int, default=600)
    parser.add_argument('--bins', type=int, default=240)
    parser.add_argument('--seed', type=int, default=20261106)
    args = parser.parse_args()
    self_check()
    if args.replicates < 1 or args.draws < 1 or args.cells < 60 or args.bins < 40 or args.bins % 4:
        parser.error('need replicates/draws >= 1, cells >= 60 and bins >= 40 divisible by 4')
    args.out.mkdir(parents=True, exist_ok=False)
    source_files = [Path(__file__).resolve()]
    import haplodose
    source_files.extend(sorted(Path(haplodose.__file__).resolve().parent.rglob('*.py')))
    provenance = dict(command=sys.argv, executable=sys.executable, python=sys.version, host=platform.node(),
                      generator_version=GENERATOR_VERSION, nuisance_model=NUISANCE_MODEL,
                      versions={name: metadata.version(name) for name in ('numpy', 'scipy', 'scikit-learn', 'jax')},
                      source_sha256={str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_files},
                      design=dict(replicates=args.replicates, draws=args.draws, cells=args.cells, bins=args.bins,
                                  seed=args.seed, condition=args.condition, null_family=6, power_family=3,
                                  null_upper_bound_alpha=.05 / 6, power_lower_bound_alpha=.05 / 3),
                      provenance_scope='runner and installed Python source; launcher must inventory native binaries/environment')
    (args.out / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    params = CloneParams(n_programmes=0, spatial_neighbors=0, seed=args.seed)
    reset = np.zeros(args.bins, bool)
    reset[::args.bins // 4] = True
    records = []
    conditions = CONDITIONS if args.condition == 'all' else (args.condition,)
    for condition in conditions:
        for r in range(args.replicates):
            rng = np.random.default_rng([args.seed, r, CONDITIONS.index(condition)])
            Y, base, phi, truth, latent = generate(rng, condition, args.cells, args.bins, return_latent=True)
            input_hashes = {name: array_hash(value) for name, value in
                            dict(counts=Y, reference=base, phi=phi, truth=truth, chromosome_reset=reset,
                                 generating_mean=latent['mean'], state_effect=latent['state_effect'],
                                 generating_depth=latent['depth']).items()}
            with (args.out / 'inputs.jsonl').open('a') as stream:
                stream.write(json.dumps(dict(condition=condition, replicate=r, input_sha256=input_hashes,
                                             generator_seed=[args.seed, r, CONDITIONS.index(condition)])) + '\n')
            settings = replace(params, seed=args.seed + r)
            observed = fit(Y, base, phi, reset, settings)
            selected = fit(Y, base, phi, reset, replace(settings, split_null_draws=args.draws))
            # Fit a single pooled parent without using its selected child groups.
            null_model = SpecimenModel(Rows(Y, np.arange(len(Y))), Y.sum(axis=1, dtype=np.float64),
                                       base[:, None], reset, phi, settings)
            pooled = null_model.profiles(Y.sum(axis=0, dtype=np.float64), null_model.L.sum(),
                                         depth_moments(null_model.L, null_model.x).sum(axis=0))
            mean = null_model.null_means(0, np.arange(len(Y)), [pooled['dose'][0]])
            obs = statistic(observed, settings.min_cells)
            null_stats = []
            for b in range(args.draws):
                sim = rng.negative_binomial(1 / phi, 1 / (1 + phi * mean)).astype(np.float32)
                replay = fit(sim, base, phi, reset, settings)
                null_stats.append(statistic(replay, settings.min_cells))
            p = (1 + sum(v >= obs for v in null_stats)) / (args.draws + 1)
            row = dict(condition=condition, replicate=r, truth_clones=int(truth.max() + 1),
                       generator_version=GENERATOR_VERSION,
                       input_sha256=input_hashes, generator_seed=[args.seed, r, CONDITIONS.index(condition)],
                       raw_clones=obs, selected_filter_clones=statistic(selected, settings.min_cells),
                       raw_membership_ari=float(adjusted_rand_score(truth, observed['R'].argmax(axis=1))),
                       selected_filter_membership_ari=float(adjusted_rand_score(truth, selected['R'].argmax(axis=1))),
                       full_replay_rank_p=p, full_replay_detected=bool(obs > 1 and p <= .05),
                       full_replay_null_clone_counts=null_stats)
            records.append(row)
            with (args.out / 'replicates.jsonl').open('a') as stream:
                stream.write(json.dumps(row) + '\n')
    summary = []
    for condition in conditions:
        rows = [r for r in records if r['condition'] == condition]
        for method in METHODS:
            detections = sum(r['full_replay_detected'] if method == 'full_replay' else r[method + '_clones'] > 1
                             for r in rows)
            lower, upper = bounds(detections, len(rows))
            family_size = 3 if condition == 'two_clones' else 6
            adjusted_lower, adjusted_upper = bounds(detections, len(rows), .05 / family_size)
            summary.append(dict(condition=condition, method=method, detections=detections, n=len(rows),
                                fraction=detections / len(rows), one_sided_95_lower=lower, one_sided_95_upper=upper,
                                simultaneous_95_lower=adjusted_lower, simultaneous_95_upper=adjusted_upper,
                                confidence_family='three alternative-condition methods' if condition == 'two_clones'
                                    else 'two null conditions x three methods',
                                confidence_rule='Bonferroni one-sided exact binomial; lower and upper are separate families',
                                endpoint='multigroup_detection_power_not_recovery' if condition == 'two_clones'
                                    else 'extra_RNA_group_rate_in_one_DNA_clone_simulations'))
    output = dict(summary=summary, params=asdict(params), seed=args.seed, draws=args.draws,
                  generator_version=GENERATOR_VERSION, nuisance_model=NUISANCE_MODEL,
                  python=sys.version, numpy=np.__version__, host=platform.node(),
                  validation_status='diagnostic_only; no pass or calibrated-clone claim',
                  design_precision=dict(null_family_size=6, power_family_size=3,
                    minimum_null_replicates_if_zero_errors=94, achieved_zero_error_resolution=args.replicates >= 94,
                    note='94 concerns simultaneous simulated-null upper bounds only, not biology or power'),
                  limitations=['simulation units, not independent biological specimens',
                               'known reference; no gene veto, alleles or spatial graph',
                               'full replay is a fitted one-parent bootstrap, not exact exchangeability',
                               'clone-count detection and descriptive membership ARI are separate endpoints'])
    (args.out / 'summary.json').write_text(json.dumps(output, indent=2) + '\n')


if __name__ == '__main__':
    main()
