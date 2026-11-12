"""Known-normal common-profile experiment; no query truth enters fitting."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.special import softmax
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def centre(x):
    return x-np.median(x)


def fingerprint(*values):
    h = hashlib.sha256()
    for value in values:
        x = np.ascontiguousarray(value)
        h.update(str((x.shape, x.dtype.str)).encode())
        h.update(x.tobytes())
    return h.hexdigest()


def estimate(Y, L, base, nuisance, phi, helper):
    """No condition or truth argument; fixed conditional NB and pooled estimates."""
    data = dict(Y=Y.astype(float), L=L, base=base, nuisance=nuisance, phi=phi,
        R=np.ones((len(Y), 1)), doses=np.zeros((1, len(base))), offset=np.zeros(len(base)))
    profile, diagnostics = helper.optimize(data)
    expected = np.sum(L[:, None]*softmax(np.log(base)+nuisance, axis=1), axis=0)
    pooled = centre(np.log((Y.sum(axis=0)+.5)/(expected+.5)))
    return profile, pooled, diagnostics


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    helper_path = Path(__file__).with_name('offset_nb_objective.py')
    helper = load('conditional_NB', helper_path)
    rec = start_run(a.out, dict(script=Path(__file__), helper=helper_path,
        design=a.task/'ANCHOR_DESIGN.md'), dict(seed_keys=[20261109, 903, 'replicate']),
        dict(scope='Known-exposure/nuisance NB experiment; independent normal calibration, biological n=0'),
        mode='reference_offset_panel')
    try:
        observations, targets, fitted, records, optimization = {}, {}, {}, [], []
        for rep in range(16):
            rng = np.random.default_rng(np.random.SeedSequence([20261109, 903, rep]))
            base = rng.lognormal(0, .5, 40)
            base /= base.sum()
            phi = np.full(40, .2)
            gamma = .25*np.sin(np.arange(40)/4)
            Ln = rng.lognormal(np.log(1000), .4, 60)
            Lq = rng.lognormal(np.log(2000), .4, 240)
            nn = (np.log(Ln)-np.log(1000))[:, None]*gamma
            nq = (np.log(Lq)-np.log(1000))[:, None]*gamma
            region = np.zeros(40)
            region[12:22] = 1.
            d = .6*np.log(1.5)*region
            technical = .35*region

            def draw(L, nuisance, shift):
                m = L[:, None]*softmax(np.log(base)+nuisance+shift, axis=1)
                return rng.negative_binomial(1/phi, 1/(1+phi*m)).astype(np.int64)

            N0, Ns = draw(Ln, nn, np.zeros(40)), draw(Ln, nn, technical)
            Q0, Qd, Qs = draw(Lq, nq, np.zeros(40)), draw(Lq, nq, d), draw(Lq, nq, d+technical)
            conditions = dict(null=(N0, Q0, np.zeros(40), np.zeros(40), np.zeros(40)),
                gain=(N0, Qd, d, np.zeros(40), np.zeros(40)),
                shared=(Ns, Qs, d, technical, technical),
                mismatched=(Ns, Qd, d, np.zeros(40), technical),
                state=(N0, Qd, np.zeros(40), d, np.zeros(40)))
            for condition, (N, Q, dna, query_other, normal_other) in conditions.items():
                key = f'{condition}__{rep}'
                nfit, npool, ndiag = estimate(N, Ln, base, nn, phi, helper)
                qfit, qpool, qdiag = estimate(Q, Lq, base, nq, phi, helper)
                profiles = dict(NB_zero=qfit, NB_anchor=centre(qfit-nfit),
                    pooled_zero=qpool, pooled_anchor=centre(qpool-npool))
                total_query = dna+query_other
                rna_zero, rna_anchor = centre(total_query), centre(total_query-normal_other)
                observations.update({key+'__'+name: value for name, value in dict(
                    N=N, Q=Q, Ln=Ln, Lq=Lq, base=base, nn=nn, nq=nq, phi=phi).items()})
                targets.update({key+'__'+name: value for name, value in dict(
                    DNA_linked_RNA=centre(dna), RNA_zero=rna_zero, RNA_anchor=rna_anchor).items()})
                optimization.append(dict(key=key, normal=ndiag, query=qdiag))
                digest = fingerprint(N, Q, Ln, Lq, base, nn, nq, phi)
                for method, profile in profiles.items():
                    fitted[key+'__'+method] = profile
                    rna = rna_anchor if method.endswith('anchor') else rna_zero
                    records.append(dict(condition=condition, replicate=rep, method=method, input_sha256=digest,
                        DNA_linked_RNA_MSE=float(np.mean((profile-centre(dna))**2)),
                        own_relative_RNA_MSE=float(np.mean((profile-rna)**2)),
                        region_contrast=float(profile[12:22].mean()-np.r_[profile[:12], profile[22:]].mean()),
                        DNA_linked_region_contrast=float(dna[12:22].mean()-np.r_[dna[:12], dna[22:]].mean())))
            for method in ('NB_zero', 'NB_anchor', 'pooled_zero', 'pooled_anchor'):
                assert np.array_equal(fitted[f'gain__{rep}__{method}'], fitted[f'state__{rep}__{method}'])
            print(json.dumps(dict(completed=rep+1, total=16)), flush=True)
        np.savez_compressed(a.out/'observations.npz', **observations)
        np.savez_compressed(a.out/'targets.npz', **targets)
        np.savez_compressed(a.out/'fitted.npz', **fitted)
        (a.out/'records.json').write_text(json.dumps(records, indent=2)+'\n')
        (a.out/'optimization.json').write_text(json.dumps(optimization, indent=2)+'\n')
        summary = []
        for condition in conditions:
            for method in profiles:
                rows = [v for v in records if v['condition'] == condition and v['method'] == method]
                summary.append(dict(condition=condition, method=method, simulation_draws=16,
                    mean_DNA_linked_RNA_MSE=float(np.mean([v['DNA_linked_RNA_MSE'] for v in rows])),
                    mean_own_relative_RNA_MSE=float(np.mean([v['own_relative_RNA_MSE'] for v in rows])),
                    mean_region_contrast=float(np.mean([v['region_contrast'] for v in rows])),
                    positive_region_contrast=sum(v['region_contrast'] > 0 for v in rows)))
        (a.out/'report.json').write_text(json.dumps(dict(summary=summary,
            optimizer_failures=sum(not d['success'] for row in optimization for d in (row['normal'], row['query'])),
            biological_specimens=0, production_changed=False, scientific_criteria_changed=False,
            scope='Known nominal exposures/nuisance; paired dependent conditions; no clone or error-control guarantee'), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU known-exposure raw NB normal/query counts and conditional profiles')


if __name__ == '__main__':
    main()
