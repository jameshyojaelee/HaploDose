"""Independent count replay, pooled arithmetic, NB gradients and recovery criteria."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import statistics

import numpy as np
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def centred(values):
    return np.asarray(values)-statistics.median(values)


def probability(base, nuisance, shift):
    logits = np.log(base)+nuisance+shift
    rates = np.exp(logits-logits.max(axis=1, keepdims=True))
    return rates/rates.sum(axis=1, keepdims=True)


def pooled(Y, L, base, nuisance):
    # Independent scalar normalization/exposure reduction and exact integer pools.
    expected = np.zeros(len(base))
    for depth, row in zip(L, nuisance):
        logits = [math.log(float(b))+float(x) for b, x in zip(base, row)]
        rates = [math.exp(v-max(logits)) for v in logits]
        total = math.fsum(rates)
        expected += [float(depth)*v/total for v in rates]
    raw = [math.log((sum(int(v) for v in Y[:, j])+.5)/(expected[j]+.5)) for j in range(len(base))]
    return centred(raw)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    helper = load('objective', root/'offset_nb_objective.py')
    reference = load('NB_reference', root/'offset_nb_numerics.py')
    original_auditor = a.task.parent/'offset_objective_20261109/stable_source_v1/benchmarks/hardening/offset_objective_audit.py'
    recovery = load('frozen_recovery_audit', original_auditor)
    folder = a.task/'anchor_v1'
    inputs = dict(script=Path(__file__), design=a.task/'ANCHOR_DESIGN.md', recovery_auditor=original_auditor,
        observations=folder/'observations.npz', targets=folder/'targets.npz', fitted=folder/'fitted.npz',
        records=folder/'records.json', summary=folder/'report.json', optimization=folder/'optimization.json',
        zero_recovery=a.task/'zero_v1/panel/summary.json')
    rec = start_run(a.out, inputs, {}, dict(scope='Known-exposure conditional counts; independent numerical and recovery audits'),
        mode='reference_offset_audit')
    try:
        _, outcomes = recovery.recovery(a.task, 'zero', 'v1')
        trace = json.loads((a.task/'zero_v1/offset_trace.json').read_text())
        assert trace and all(v['max_previous_offset'] == 0. for v in trace)
        observed, target, fitted = (np.load(inputs[k]) for k in ('observations', 'targets', 'fitted'))
        records = json.loads(inputs['records'].read_text())
        report = json.loads(inputs['summary'].read_text())
        assert len(records) == 320 and len(report['summary']) == 20
        eps = np.finfo(float).eps
        metrics, gradient_checks, stationarity = [], [], []
        counts_replayed = 0
        for rep in range(16):
            rng = np.random.default_rng(np.random.SeedSequence([20261109, 903, rep]))
            base = rng.lognormal(0, .5, 40)
            base /= base.sum()
            Ln, Lq = rng.lognormal(np.log(1000), .4, 60), rng.lognormal(np.log(2000), .4, 240)
            gamma, phi = .25*np.sin(np.arange(40)/4), np.full(40, .2)
            nn, nq = (np.log(Ln)-np.log(1000))[:, None]*gamma, (np.log(Lq)-np.log(1000))[:, None]*gamma
            r = np.zeros(40)
            r[12:22] = 1.
            gain, technical = .6*np.log(1.5)*r, .35*r

            def sample(L, nuisance, shift):
                means = L[:, None]*probability(base, nuisance, shift)
                return rng.negative_binomial(1/phi, 1/(1+phi*means)).astype(np.int64)

            N0, Ns = sample(Ln, nn, np.zeros(40)), sample(Ln, nn, technical)
            Q0, Qd, Qs = sample(Lq, nq, np.zeros(40)), sample(Lq, nq, gain), sample(Lq, nq, gain+technical)
            cases = dict(null=(N0, Q0, np.zeros(40), np.zeros(40), np.zeros(40)),
                gain=(N0, Qd, gain, np.zeros(40), np.zeros(40)),
                shared=(Ns, Qs, gain, technical, technical),
                mismatched=(Ns, Qd, gain, np.zeros(40), technical),
                state=(N0, Qd, np.zeros(40), gain, np.zeros(40)))
            for condition, (N, Q, dna, qother, nother) in cases.items():
                key = f'{condition}__{rep}'
                for name, value in dict(N=N, Q=Q, Ln=Ln, Lq=Lq, base=base, nn=nn, nq=nq, phi=phi).items():
                    assert np.array_equal(observed[key+'__'+name], value), (key, name)
                counts_replayed += N.size+Q.size
                digest = hashlib.sha256()
                for value in (N, Q, Ln, Lq, base, nn, nq, phi):
                    value = np.ascontiguousarray(value)
                    digest.update(str((value.shape, value.dtype.str)).encode())
                    digest.update(value.tobytes())
                expected_targets = dict(DNA_linked_RNA=centred(dna), RNA_zero=centred(dna+qother),
                    RNA_anchor=centred(dna+qother-nother))
                for name, value in expected_targets.items():
                    assert np.array_equal(target[key+'__'+name], value)
                independent_pooled = dict(pooled_zero=pooled(Q, Lq, base, nq))
                independent_pooled['pooled_anchor'] = centred(independent_pooled['pooled_zero']-pooled(N, Ln, base, nn))
                for method, value in independent_pooled.items():
                    budget = 128*(N.size+Q.size)*eps
                    assert np.max(np.abs(fitted[key+'__'+method]-value)) <= budget
                for method in ('NB_zero', 'NB_anchor', 'pooled_zero', 'pooled_anchor'):
                    estimate = fitted[key+'__'+method]
                    dna_mse = math.fsum((float(x)-float(y))**2 for x, y in zip(estimate, centred(dna)))/40
                    rna = expected_targets['RNA_anchor' if method.endswith('anchor') else 'RNA_zero']
                    rna_mse = math.fsum((float(x)-float(y))**2 for x, y in zip(estimate, rna))/40
                    contrast = math.fsum(estimate[12:22])/10-math.fsum(np.r_[estimate[:12], estimate[22:]])/30
                    row = next(v for v in records if v['condition'] == condition and v['replicate'] == rep and v['method'] == method)
                    assert row['input_sha256'] == digest.hexdigest()
                    budget = 128*40*eps*max(1., dna_mse, rna_mse, abs(contrast))
                    assert abs(row['DNA_linked_RNA_MSE']-dna_mse) <= budget
                    assert abs(row['own_relative_RNA_MSE']-rna_mse) <= budget
                    assert abs(row['region_contrast']-contrast) <= budget
                    metrics.append(row)
                qfit = fitted[key+'__NB_zero']
                nfit = centred(qfit-fitted[key+'__NB_anchor'])
                for role, Y, L, nuisance, estimate in (('normal', N, Ln, nn, nfit), ('query', Q, Lq, nq, qfit)):
                    data = dict(Y=Y.astype(float), L=L, base=base, nuisance=nuisance, phi=phi,
                        R=np.ones((len(Y), 1)), doses=np.zeros((1, 40)), offset=np.zeros(40))
                    value, gradient, budget = helper.gain_gradient(data, estimate)
                    assert value >= -budget
                    stationarity.append(dict(key=key, role=role,
                        normalized_gradient=float(np.max(np.abs(gradient))/Y.sum())))
                    if rep == 0:
                        trial = estimate.copy()
                        trial[[0, 12, 39]] += [.1, -.1, .05]
                        analytic = helper.gain_gradient(data, trial)[1]
                        for j in (0, 12, 39):
                            h = eps**(1/3)
                            differences = []
                            for step in (h, h/2, h/4):
                                plus, minus = trial.copy(), trial.copy()
                                plus[j] += step
                                minus[j] -= step
                                differences.append((reference.full_q(data, plus)-reference.full_q(data, minus))/(2*step))
                            gradient_budget = 8*max(abs(differences[1]-differences[0]), abs(differences[2]-differences[1]))
                            gradient_budget += 128*eps*max(1., abs(reference.full_q(data, trial)))/(h/4)
                            assert abs(analytic[j]-differences[-1]) <= gradient_budget
                            gradient_checks.append(dict(key=key, role=role, bin=j, gradient=float(analytic[j]),
                                independent=float(differences[-1]), budget=gradient_budget))
            for method in ('NB_zero', 'NB_anchor', 'pooled_zero', 'pooled_anchor'):
                assert np.array_equal(fitted[f'gain__{rep}__{method}'], fitted[f'state__{rep}__{method}'])
            for name in ('N', 'Q', 'Ln', 'Lq', 'base', 'nn', 'nq', 'phi'):
                assert np.array_equal(observed[f'gain__{rep}__{name}'], observed[f'state__{rep}__{name}'])
        for row in report['summary']:
            rows = [v for v in metrics if v['condition'] == row['condition'] and v['method'] == row['method']]
            for name in ('DNA_linked_RNA_MSE', 'own_relative_RNA_MSE', 'region_contrast'):
                expected = math.fsum(v[name] for v in rows)/16
                assert abs(row['mean_'+name]-expected) <= 128*eps*max(1., abs(expected))
            assert row['positive_region_contrast'] == sum(v['region_contrast'] > 0 for v in rows)
        (a.out/'report.json').write_text(json.dumps(dict(all_checks_passed=True,
            zero_recovery=outcomes, independent_recovery_criteria=30, zero_offset_updates=len(trace),
            exact_count_entries_replayed=counts_replayed, case_method_metrics=320,
            gradient_checks=gradient_checks, stationarity=stationarity, summary=report['summary'],
            biological_specimens=0, production_changed=False, scientific_criteria_changed=False,
            scope='Known exposure/nuisance model; no biological, clone-error or uncertainty guarantee'), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU independent sampler/pooled arithmetic/NB finite differences/recovery recount')


if __name__ == '__main__':
    main()
