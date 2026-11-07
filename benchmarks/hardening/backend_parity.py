"""Bounded installed-wheel CPU/GPU checks with independent numerical references.

Run on a compute node. Each device runs in a separate process; GPU failure never
falls back to CPU. Discrete boundary changes are reported separately from numeric
disagreement. This is a computational diagnostic, not biological validation.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys
import zipfile

import numpy as np
from scipy import integrate, optimize, special, stats

from haplodose.io import new_output_dir, write_json
from haplodose.provenance import file_record, finish_run, start_run


def fresh_output(task_root, output):
    root, output = task_root.resolve(strict=True), output.resolve()
    relative = output.relative_to(root)
    if not relative.parts or relative.parts[0] == 'baseline':
        raise ValueError('Output must be a new task_root descendant outside baseline')
    if output.exists():
        raise FileExistsError(output)
    return new_output_dir(output)


def verify_installed_wheel(wheel):
    import haplodose

    package = Path(haplodose.__file__).resolve().parent
    checkout = Path(__file__).resolve().parents[2] / 'python'
    if package.is_relative_to(checkout):
        raise RuntimeError('Backend diagnostic requires an installed wheel, not the source checkout')
    checked = {}
    with zipfile.ZipFile(wheel) as archive:
        for name in archive.namelist():
            member = Path(name)
            if member.parts[0] != 'haplodose' or member.suffix not in ('.py', '.so', '.tsv'):
                continue
            installed = package.joinpath(*member.parts[1:])
            expected = hashlib.sha256(archive.read(name)).hexdigest()
            actual = file_record(installed)['sha256']
            if actual != expected:
                raise RuntimeError(f'Installed package differs from the supplied wheel: {name}')
            checked[name] = actual
    if not checked:
        raise RuntimeError('Supplied wheel has no HaploDose package sources')
    return dict(package_path=str(package), package_version=haplodose.__version__, verified_files=checked)


def forward_backward64(emit, trans, init, reset):
    """Independent NumPy/SciPy log recursion with chromosome resets."""
    T, S = emit.shape
    alpha, beta = np.empty((T, S)), np.zeros((T, S))
    total = 0.
    for t in range(T):
        prior = init if reset[t] else special.logsumexp(alpha[t - 1, :, None] + trans, axis=0)
        value = prior + emit[t]
        scale = special.logsumexp(value)
        alpha[t] = value - scale
        total += scale
    for t in range(T - 2, -1, -1):
        if not reset[t + 1]:
            value = special.logsumexp(trans + emit[t + 1][None] + beta[t + 1][None], axis=1)
            beta[t] = value - value.max()
    post = special.softmax(alpha + beta, axis=1)
    return post, float(total)


def pln_quadrature(y, log_mu, tau):
    """Scaled float64 SciPy integral, following tests/test_numerics.py's reference."""
    log_density = lambda offset: stats.poisson.logpmf(y, np.exp(log_mu + offset)) + stats.norm.logpdf(offset, 0, tau)
    mode = optimize.minimize_scalar(lambda o: -log_density(o), bounds=(-12 * tau, 12 * tau), method='bounded').x
    top = log_density(mode)
    integral, error = integrate.quad(lambda o: np.exp(log_density(o) - top), mode - 12 * tau,
                                     mode + 12 * tau, points=[mode], epsabs=0., limit=400)
    return float(top + np.log(integral)), float(error / integral)


def programme_objective(u, y, responsibility, log_rate, phi, W, precision):
    mean = y.sum() * special.softmax(log_rate + u @ W, axis=1)
    values = stats.nbinom.logpmf(y[None], 1 / phi, 1 / (1 + phi * mean)).sum(axis=1)
    return float(responsibility @ values - .5 * np.sum(precision * u ** 2))


def run_worker(args, out):
    import jax
    import jax.numpy as jnp
    from haplodose.device import select
    from haplodose.kernels import estep, hmm
    from haplodose import profile
    from haplodose.kernels.likelihood import pln_logpmf

    jax.config.update('jax_enable_x64', False)
    backend = select(args.worker)
    seen = {}

    def observe(name, function):
        def checked(*positional, **keywords):
            result = function(*positional, **keywords)
            arrays = result if isinstance(result, tuple) else (result,)
            devices = sorted({(d.platform, d.id) for array in arrays if hasattr(array, 'devices') for d in array.devices()})
            if not devices or any(platform != args.worker for platform, _ in devices):
                raise RuntimeError(f'{name} executed on {devices}, expected {args.worker}')
            for array in arrays:
                if hasattr(array, 'block_until_ready'):
                    array.block_until_ready()
            item = seen.setdefault(name, dict(calls=0, devices=devices))
            item['calls'] += 1
            return result
        return checked

    # Observation wrappers preserve the actual installed implementation and record
    # device placement at its JAX outputs, rather than relying on an import flag.
    profile.forward_backward_batch = observe('profile_forward_backward', profile.forward_backward_batch)
    profile.viterbi_batch = observe('profile_viterbi', profile.viterbi_batch)
    profile.pln_logpmf = observe('profile_pln', profile.pln_logpmf)
    estep._chunk_loglik = observe('cell_nb_loglik', estep._chunk_loglik)
    estep._chunk_scores_fit = observe('programme_fit', estep._chunk_scores_fit)
    rng = np.random.default_rng(args.seed)
    J, F, K, n = 120, 2, 3, 3
    arrays, result = {}, {}

    # PLN reference uses the actual float32 input values as well as the intended
    # float64 values, separating input quantization from quadrature/kernel error.
    pln_rows = []
    for y, mean, tau in [(0., .3, .2), (7., 3., .2), (300., 450., .2), (5e4, 4e4, .2), (1e7, 1e7, .15)]:
        log_mu = np.log(mean)
        got = observe('standalone_pln', pln_logpmf)(y, log_mu, tau)
        reference, quad_error = pln_quadrature(float(np.float32(y)), float(np.float32(log_mu)), float(np.float32(tau)))
        intended, intended_error = pln_quadrature(y, log_mu, tau)
        value = float(got)
        pln_rows.append(dict(y=y, mean=mean, tau=tau, value=value, input_matched_reference=reference,
                             intended_float64_reference=intended, quadrature_relative_error=quad_error,
                             intended_quadrature_relative_error=intended_error,
                             kernel_absolute_error=abs(value - reference),
                             input_quantization_effect=reference - intended,
                             output_ulp=float(abs(np.spacing(np.float32(value))))))
    result['pln'] = pln_rows

    prior = np.array([.2, .5, .2, .1])
    transition = .002 * prior[None] / (1 - prior[:, None])
    np.fill_diagonal(transition, 1 - .002)
    emit = rng.normal(0, .3, (J, 4))
    emit[np.arange(J), np.repeat([1, 2, 1, 0], J // 4)] += 1.
    emit, trans, init = emit.astype(np.float32), np.log(transition).astype(np.float32), np.log(prior).astype(np.float32)
    reset = np.zeros(J, bool)
    reset[::J // 4] = True
    reference_post, reference_ll = forward_backward64(emit.astype(float), trans.astype(float), init.astype(float), reset)
    result['hmm'] = {}
    arrays['hmm_reference_post'] = reference_post
    for method in ('scan', 'assoc'):
        post, ll = observe('hmm_' + method, lambda *a: hmm.forward_backward(*a, method=method))(
            jnp.asarray(emit), jnp.asarray(trans), jnp.asarray(init), jnp.asarray(reset))
        actual = np.asarray(post)
        arrays['hmm_' + method + '_post'] = actual
        result['hmm'][method] = dict(max_posterior_error=float(np.max(np.abs(actual - reference_post))),
                                    log_likelihood=float(ll), reference_log_likelihood=reference_ll,
                                    posterior_call_flips_at_095=int(np.sum((actual >= .95) != (reference_post >= .95))),
                                    map_flips_vs_reference=int(np.sum(actual.argmax(axis=1) != reference_post.argmax(axis=1))),
                                    smallest_reference_map_margin=float(np.min(np.diff(np.sort(reference_post, axis=1)[:, -2:], axis=1))))

    short_emit = rng.normal(0, 2, (7, 3)).astype(np.float32)
    short_trans = np.log(rng.dirichlet(np.ones(3), 3)).astype(np.float32)
    short_init = np.log(rng.dirichlet(np.ones(3))).astype(np.float32)
    short_reset = np.array([True, False, False, True, False, False, False])
    paths = np.asarray(list(itertools.product(range(3), repeat=7)))
    scores = np.array([sum(float(short_emit[t, path[t]]) + float(short_init[path[t]] if short_reset[t]
        else short_trans[path[t - 1], path[t]]) for t in range(7)) for path in paths])
    best = np.argsort(scores)[-2:]
    short_path = np.asarray(observe('standalone_viterbi', hmm.viterbi)(jnp.asarray(short_emit),
        jnp.asarray(short_trans), jnp.asarray(short_init), jnp.asarray(short_reset)))
    # Path-score reduction bound depends on this chain's length and magnitudes.
    eps = np.finfo(np.float32).eps
    score_bound = 14 * eps / (1 - 14 * eps) * 7 * (
        np.max(np.abs(short_emit)) + max(np.max(np.abs(short_trans)), np.max(np.abs(short_init))))
    result['hmm']['enumeration'] = dict(path_disagreements=int(np.sum(short_path != paths[best[-1]])),
        best_vs_second_score_margin=float(scores[best[-1]] - scores[best[-2]]),
        arithmetic_path_score_budget=float(score_bound),
        best_path_separated_beyond_roundoff=bool(scores[best[-1]] - scores[best[-2]] > 2 * score_bound))

    log_rate = rng.normal(-np.log(J), .3, (K, J))
    phi = rng.uniform(.05, .5, J)
    W, precision = rng.normal(0, .25, (F, J)), np.array([.8, 1.1])
    responsibility = np.array([[.15, .55, .30], [.5, .25, .25], [.05, .1, .85]])
    y = np.vstack([rng.negative_binomial(1 / phi, 1 / (1 + phi * 2000 * special.softmax(log_rate[k])) ) for k in range(n)])
    L = y.sum(axis=1)
    likelihood = estep.cell_loglik(y, L, log_rate, phi)
    mean = L[:, None, None] * special.softmax(log_rate, axis=1)[None]
    ref_likelihood = stats.nbinom.logpmf(y[:, None], 1 / phi, 1 / (1 + phi * mean)).sum(axis=2)
    ref_likelihood -= ref_likelihood[:, :1]
    encoded_log_rate = log_rate[:1].astype(np.float32).astype(float) + (log_rate - log_rate[:1]).astype(np.float32).astype(float)
    encoded_size = (1 / phi).astype(np.float32).astype(float)
    encoded_mean = L[:, None, None] * special.softmax(encoded_log_rate, axis=1)[None]
    encoded_reference = stats.nbinom.logpmf(y[:, None], encoded_size,
        encoded_size / (encoded_size + encoded_mean)).sum(axis=2)
    encoded_reference -= encoded_reference[:, :1]
    arrays['cell_loglik'], arrays['cell_loglik_reference'] = likelihood, ref_likelihood
    arrays['cell_loglik_input_matched_reference'] = encoded_reference
    result['cell_nb'] = dict(max_reference_error=float(np.max(np.abs(likelihood - ref_likelihood))),
                             max_input_matched_reference_error=float(np.max(np.abs(likelihood - encoded_reference))),
                             max_input_quantization_effect=float(np.max(np.abs(encoded_reference - ref_likelihood))),
                             map_flips_vs_reference=int(np.sum(likelihood.argmax(axis=1) != ref_likelihood.argmax(axis=1))),
                             reference_assignment_margins=np.diff(np.sort(ref_likelihood, axis=1)[:, -2:], axis=1).ravel().tolist())
    u, diagnostics = estep.programme_scores(y, L, responsibility, log_rate, phi, W, precision, return_diagnostics=True)
    arrays['programme_scores'] = u
    result['programme'] = dict(diagnostics={k: v.tolist() for k, v in diagnostics.items()}, reference=[])
    for c in range(n):
        objective = lambda score: programme_objective(score, y[c], responsibility[c].astype(np.float32).astype(float),
            log_rate.astype(np.float32).astype(float), 1 / encoded_size, W.astype(np.float32).astype(float),
            precision.astype(np.float32).astype(float))
        candidates = [optimize.minimize(lambda score: -objective(score), start, method='BFGS')
                      for start in (np.zeros(F), u[c].astype(float), -u[c].astype(float))]
        best = min(candidates, key=lambda fitted: fitted.fun)
        result['programme']['reference'].append(dict(cell=c, local_optimizer_success=bool(best.success),
            local_optimizer_message=str(best.message), local_gradient_norm=float(np.linalg.norm(best.jac)),
            objective_gain_from_zero=float(objective(u[c]) - objective(np.zeros(F))),
            independent_local_objective_gap=float(objective(best.x) - objective(u[c])),
            max_score_disagreement=float(np.max(np.abs(best.x - u[c]))),
            reference_scope='float32-input-matched working likelihood; best of three local BFGS starts, not global proof'))

    expected = np.full(J, 5000.)
    injected = np.repeat([0., .6 * np.log(1.5), 0., .6 * np.log(.5)], J // 4)
    counts = rng.poisson(expected * np.exp(injected + rng.normal(0, .08, J))).astype(float)
    variance = np.linspace(0, .05 ** 2, J)
    params = profile.ProfileParams()
    center, beta, tau = profile.fit_pooled(counts, expected, reset, params, v=variance)
    post, path, dose, _ = profile.segment(counts[None], expected[None], reset, [center], beta,
                                         np.sqrt(tau ** 2 + variance)[None], params)
    segments = profile.segments_table(post[0], path[0], dose[0], counts, expected, reset, beta, tau, params,
                                      reference_var=variance)
    arrays.update(profile_post=post[0], profile_path=path[0], profile_dose=dose[0])
    arrays['profile_segment_se_per_bin'] = np.empty(J)
    arrays['profile_resolved_per_bin'] = np.empty(J, np.int32)
    for row in segments:
        sl = slice(row['first_bin'], row['last_bin'] + 1)
        arrays['profile_segment_se_per_bin'][sl] = row['se_rna_log_ratio_approx']
        arrays['profile_resolved_per_bin'][sl] = row['resolved']
    result['profile'] = dict(center=center, beta=beta, tau=tau, segments=segments,
                              uncertainty_scope='independent-bin conditional approximation, uncalibrated',
                              default_hmm_method=hmm.default_method())
    if not all(np.isfinite(a).all() for a in arrays.values()):
        raise FloatingPointError('Nonfinite inference/reference output')
    result.update(requested_backend=args.worker, actual_backend=backend, inference_devices=seen,
                  jax_version=jax.__version__, jax_enable_x64=bool(jax.config.jax_enable_x64),
                  visible_devices=[str(d) for d in jax.devices()], array_dtypes={k: str(v.dtype) for k, v in arrays.items()})
    import resource
    result['cpu_max_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    result['device_memory_stats'] = {}
    for device in jax.devices():
        try:
            result['device_memory_stats'][str(device)] = device.memory_stats()
        except (AttributeError, RuntimeError, NotImplementedError) as error:
            result['device_memory_stats'][str(device)] = {'unavailable': str(error)}
    np.savez_compressed(out / 'arrays.npz', **arrays)
    write_json(result, out / 'summary.json')
    return result


def comparison(out, devices):
    if set(devices) != {'cpu', 'gpu'}:
        return dict(status='single_backend_reference_checks_only', cross_backend_parity='not_run')
    a, b = np.load(out / 'cpu' / 'arrays.npz'), np.load(out / 'gpu' / 'arrays.npz')
    result = {}
    for key in a.files:
        if a[key].shape != b[key].shape:
            result[key] = dict(shape_agrees=False, cpu_shape=list(a[key].shape), gpu_shape=list(b[key].shape))
        else:
            result[key] = dict(shape_agrees=True, max_absolute_difference=float(np.max(np.abs(a[key] - b[key]))),
                               exact_equality=bool(np.array_equal(a[key], b[key])))
    result['decision_changes'] = dict(
        profile_viterbi_bins=int(np.sum(a['profile_path'] != b['profile_path'])),
        profile_posterior_calls_at_095=int(np.sum((a['profile_post'] >= .95) != (b['profile_post'] >= .95))),
        profile_resolved_call_bins=int(np.sum(a['profile_resolved_per_bin'] != b['profile_resolved_per_bin'])),
        cell_map_assignments=int(np.sum(a['cell_loglik'].argmax(axis=1) != b['cell_loglik'].argmax(axis=1))))
    result['assessment'] = 'Report independent-reference errors and discrete changes; no universal epsilon or biological stability claim'
    a.close()
    b.close()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--wheel', type=Path, required=True)
    parser.add_argument('--seed', type=int, default=20261106)
    parser.add_argument('--devices', nargs='+', choices=['cpu', 'gpu'], default=['cpu', 'gpu'])
    parser.add_argument('--worker', choices=['cpu', 'gpu'], help=argparse.SUPPRESS)
    args = parser.parse_args()
    installed = verify_installed_wheel(args.wheel)
    out = fresh_output(args.task_root, args.out)
    record = start_run(out, dict(probe_source=__file__, wheel=args.wheel),
                       dict(seed=args.seed, devices=args.devices, worker=args.worker, installed=installed),
                       dict(scope='synthetic computational backend checks', biological_validation='not_established'),
                       mode='backend_parity_worker' if args.worker else 'backend_parity')
    try:
        if args.worker:
            summary = run_worker(args, out)
            finish_run(out, record, backend=summary['actual_backend'])
        else:
            statuses = {}
            for device in dict.fromkeys(args.devices):
                command = [sys.executable, str(Path(__file__).resolve()), '--task-root', str(args.task_root.resolve()),
                           '--out', str(out / device), '--wheel', str(args.wheel.resolve()), '--seed', str(args.seed),
                           '--worker', device]
                environment = dict(os.environ)
                environment.setdefault('XLA_PYTHON_CLIENT_PREALLOCATE', 'false')
                with (out / (device + '.stdout.txt')).open('w') as stdout, (out / (device + '.stderr.txt')).open('w') as stderr:
                    completed = subprocess.run(command, env=environment, stdout=stdout, stderr=stderr)
                statuses[device] = completed.returncode
            write_json(dict(worker_exit_codes=statuses,
                            parity=comparison(out, list(statuses)) if not any(statuses.values()) else None),
                       out / 'comparison.json')
            if any(statuses.values()):
                raise RuntimeError(f'Backend workers failed: {statuses}; inspect separate stderr and manifests')
            record['children'] = {device: json.loads((out / device / 'run_manifest.json').read_text())['output']['id']
                                  for device in statuses}
            finish_run(out, record, backend='separate_processes')
        print(json.dumps({'out': str(out), 'status': 'complete'}))
    except Exception as error:
        finish_run(out, record, error=error)
        raise


if __name__ == '__main__':
    main()
