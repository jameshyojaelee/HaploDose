"""Unpromoted conditional offset candidate; its complete workflow fails recovery.

Exact conditional-pool arithmetic does not make this a compatible replacement
for the package's legacy clone-only shared-offset/profile procedure.
"""
import numpy as np
from scipy.special import logsumexp


def update_offset(model, R, doses, keep):
    Rk = np.asarray(R[:, keep], np.float64)
    lr = np.log(model.base)[None]+np.asarray(doses)
    observed, expected = np.zeros(model.J), np.zeros(model.J)
    for a, b, Yc in model.rows.chunks():
        weights = Rk[a:b]
        observed += (weights.sum(axis=1)[:, None]*Yc).sum(axis=0)
        nuisance = model.x[a:b, None]*model.gamma[None]+model.U[a:b]@model.W
        for k in range(len(doses)):
            eta = lr[k]+nuisance
            p = np.exp(eta-logsumexp(eta, axis=1, keepdims=True))
            expected += (weights[:, k]*model.L[a:b])@p
    residual = np.log((observed+.5)/(expected+.5))
    model.offset = residual-np.median(residual)
