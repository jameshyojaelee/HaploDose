"""Research-only latent-weight dispersion moments; production remains untouched."""
from contextlib import contextmanager

import numpy as np
from haplodose.clones import SpecimenModel


@contextmanager
def weighted_dispersion(trace, save_moments=True):
    original = SpecimenModel._update_phi

    def update(model, R, doses):
        R = np.asarray(R, float)
        if R.shape != (model.rows.n, len(doses)) or not np.isfinite(R).all() or (R < 0).any():
            raise ValueError('Aligned finite nonnegative latent weights required')
        lr = model.log_rates(doses)
        rates = np.exp(lr)
        num, den, expanded, absolute = [np.zeros(model.J) for _ in range(4)]
        for a, b, Yc in model.rows.chunks():
            shift = np.exp(model.x[a:b, None]*model.gamma[None]+model.U[a:b]@model.W)
            first, second = np.zeros(Yc.shape), np.zeros(Yc.shape)
            mass = R[a:b].sum(axis=1)[:, None]
            for k in range(len(doses)):
                rc = rates[k]*shift
                m = model.L[a:b, None]*rc/rc.sum(axis=1, keepdims=True)
                w = R[a:b, k, None]
                num += (w*((Yc-m)**2-m)).sum(axis=0)
                den += (w*m*m).sum(axis=0)
                first += w*m
                second += w*m*m
            expanded += (mass*Yc.astype(float)**2-2*Yc*first+second-first).sum(axis=0)
            absolute += (mass*Yc.astype(float)**2+2*Yc*first+second+first).sum(axis=0)
        operations = 128*model.rows.n*len(doses)
        unit = np.finfo(float).eps
        budget = operations*unit/(1-operations*unit)
        assert (np.abs(num-expanded) <= budget*np.maximum(1., absolute)).all()
        phi = np.clip(np.divide(num, den, out=np.full(model.J, .05), where=den > 0), 1e-3, 10.)
        result = .5*phi+.5*np.median(phi)
        before = model.phi.copy()
        original(model, R, doses)
        counterfactual = model.phi.copy()
        model.phi = before
        record = dict(cells=model.rows.n, components=len(doses),
            total_latent_mass=float(R.sum()), min_row_mass=float(R.sum(axis=1).min()),
            max_row_mass=float(R.sum(axis=1).max()),
            soft_phi_median=float(np.median(result)), MAP_same_inputs_phi_median=float(np.median(counterfactual)),
            max_phi_change=float(np.max(np.abs(result-counterfactual))),
            max_moment_expansion_difference=float(np.max(np.abs(num-expanded))), numerical_factor=budget)
        if save_moments:
            record.update(numerator=num.tolist(), denominator=den.tolist(), soft_phi=result.tolist(),
                          MAP_same_inputs_phi=counterfactual.tolist())
        trace.append(record)
        model.phi = result

    SpecimenModel._update_phi = update
    try:
        yield
    finally:
        SpecimenModel._update_phi = original
