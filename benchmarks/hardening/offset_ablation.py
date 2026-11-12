"""Separate retained-population and nuisance-mean changes without promoting either."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from scipy.special import logsumexp
from haplodose.clones import SpecimenModel
from haplodose.provenance import start_run, finish_run


@contextmanager
def ablation(mode, trace):
    original = SpecimenModel._update_offset

    def update(model, R, doses, keep):
        Rk = np.asarray(R[:, keep], float)
        lr = np.log(model.base)[None]+np.asarray(doses)
        if mode == 'population_only':
            N = Rk.T@model.L
            rates = np.exp(lr-logsumexp(lr, axis=1, keepdims=True))
            expected = (N[:, None]*rates).sum(axis=0)
            observed = model.rows.weighted_sum(Rk.sum(axis=1)[:, None])[0]
        elif mode == 'nuisance_only':
            expected = np.zeros(model.J)
            for a, b, _ in model.rows.chunks():
                nuisance = model.x[a:b, None]*model.gamma[None]+model.U[a:b]@model.W
                for k in range(len(doses)):
                    eta = lr[k]+nuisance
                    p = np.exp(eta-logsumexp(eta, axis=1, keepdims=True))
                    expected += (Rk[a:b, k]*model.L[a:b])@p
            observed = model.Ytot
        else:
            raise ValueError(mode)
        o = np.log((observed+.5)/(expected+.5))
        model.offset = o-np.median(o)
        trace.append(dict(mode=mode, kept_components=len(keep),
            removed_components=R.shape[1]-len(keep), min_kept_mass=float(Rk.sum(axis=1).min()),
            max_abs_offset=float(np.max(np.abs(model.offset)))))
    SpecimenModel._update_offset = update
    try:
        yield
    finally:
        SpecimenModel._update_offset = original


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--mode', choices=['population_only', 'nuisance_only'], required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    path = Path(__file__).with_name('membership_stage_panel.py')
    spec = importlib.util.spec_from_file_location('recovery', path)
    recovery = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(recovery)
    rec = start_run(a.out, dict(script=Path(__file__), runner=path, design=a.task/'ABLATION_DESIGN.md'),
        dict(mode=a.mode), dict(scope='Failed candidate component diagnosis, no default promotion'), mode='offset_ablation')
    trace, previous = [], sys.argv
    try:
        sys.argv = [str(path), '--out', str(a.out/'panel'), '--stage', 'hard',
            '--panel-root', str(a.task.parent/'recovery_20261107/recovery_panel_v1'),
            '--assertion-source', str(a.task.parent/'recovery_20261107/baseline/tests/test_split_mode.py')]
        with ablation(a.mode, trace):
            recovery.main()
        (a.out/'offset_trace.json').write_text(json.dumps(trace, indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finally:
        sys.argv = previous
    finish_run(a.out, rec, backend='CPU fixed original recovery with one offset component ablated')


if __name__ == '__main__':
    main()
