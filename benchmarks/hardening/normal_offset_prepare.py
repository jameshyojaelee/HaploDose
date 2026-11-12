"""Freeze fresh NB observation panels separately from truth and leakage challenges."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import softmax
from haplodose.provenance import start_run, finish_run


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--seed-key', type=int, default=904)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    rec = start_run(a.out, dict(script=Path(__file__), design=a.task/'DESIGN.md'),
        dict(seed_keys=[20261109, a.seed_key, 'replicate']),
        dict(scope='Known-exposure/nuisance conditional NB panel, biological n=0'), mode='normal_offset_prepare')
    try:
        observed, targets, keys = {}, {}, []
        for rep in range(16):
            rng = np.random.default_rng(np.random.SeedSequence([20261109, a.seed_key, rep]))
            base = rng.lognormal(0, .5, 40)
            base /= base.sum()
            phi = np.full(40, .2)
            gamma = .25*np.sin(np.arange(40)/4)
            Ln, Lq = rng.lognormal(np.log(1000), .4, 60), rng.lognormal(np.log(2000), .4, 240)
            nn = (np.log(Ln)-np.log(1000))[:, None]*gamma
            nq = (np.log(Lq)-np.log(1000))[:, None]*gamma
            region = np.zeros(40)
            region[12:22] = 1.
            gain, technical = .6*np.log(1.5)*region, .35*region

            def draw(L, nuisance, shift):
                m = L[:, None]*softmax(np.log(base)+nuisance+shift, axis=1)
                return rng.negative_binomial(1/phi, 1/(1+phi*m)).astype(np.int64)

            N0, Ns = draw(Ln, nn, np.zeros(40)), draw(Ln, nn, technical)
            Q0, Qd, Qs = draw(Lq, nq, np.zeros(40)), draw(Lq, nq, gain), draw(Lq, nq, gain+technical)
            conditions = dict(null=(N0, Q0, np.zeros(40), np.zeros(40), np.zeros(40)),
                gain=(N0, Qd, gain, np.zeros(40), np.zeros(40)),
                shared=(Ns, Qs, gain, technical, technical),
                mismatched=(Ns, Qd, gain, np.zeros(40), technical),
                state=(N0, Qd, np.zeros(40), gain, np.zeros(40)))
            for condition, (N, Q, dna, qother, nother) in conditions.items():
                key = f'{condition}__{rep}'
                keys.append(key)
                observed.update({key+'__'+name: value for name, value in dict(
                    N=N, Q=Q, Ln=Ln, Lq=Lq, base=base, nn=nn, nq=nq, phi=phi).items()})
                centre = lambda x: x-np.median(x)
                targets.update({key+'__'+name: value for name, value in dict(
                    DNA_linked_RNA=centre(dna), RNA_zero=centre(dna+qother),
                    RNA_anchor=centre(dna+qother-nother)).items()})
        np.savez_compressed(a.out/'observations.npz', **observed)
        np.savez_compressed(a.out/'targets.npz', **targets)
        changed = {k: value[:, ::-1]+7 if k.endswith('__Q') else value for k, value in observed.items()}
        np.savez_compressed(a.out/'observations_query_changed.npz', **changed)
        np.savez_compressed(a.out/'targets_changed.npz', **{k: -3*v for k, v in targets.items()})
        (a.out/'metadata.json').write_text(json.dumps(dict(keys=keys, normal_fit_rows='even',
            normal_validation_rows='odd', reverse_fold_also_scored=True, biological_specimens=0,
            query_perturbation='reverse query bin order and add7; no performance interpretation',
            truth_perturbation='multiply all truth effects by-3; no performance interpretation'), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU raw NB panel and sealed independent leakage challenges')


if __name__ == '__main__':
    main()
