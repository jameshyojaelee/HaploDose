"""Scalar expected-information correction; a local isotropic approximation only."""
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import softmax


def information(L, base, nuisance, phi):
    p = softmax(np.log(base)+nuisance, axis=1)
    m = L[:, None]*p
    w = m/(1+phi*m)
    trace = w.sum(axis=1)*(1+(p*p).sum(axis=1))-2*(p*w).sum(axis=1)
    value = float(trace.sum()/(len(base)-1))
    if not np.isfinite(value) or value <= 0:
        raise ValueError('Positive finite information on the identifiable logit space is required')
    return value


def fraction(parameter, value, mean_fold_information):
    return float(value*parameter/(value*parameter+mean_fold_information*(1-parameter)))


def selection(N, L, base, nuisance, phi, helper, profiles):
    """Select prior strength with normals alone, using additive full-refit information."""
    folds = [np.arange(parity, len(N), 2) for parity in (0, 1)]
    infos = [information(L[rows], base, nuisance[rows], phi) for rows in folds]
    mean_info, full_info = sum(infos)/2, sum(infos)
    actual_full = information(L, base, nuisance, phi)
    budget = 128*len(N)*len(base)*np.finfo(float).eps*max(1., full_info)
    assert abs(full_info-actual_full) <= budget
    validation = [dict(Y=N[rows].astype(float), L=L[rows], base=base, nuisance=nuisance[rows], phi=phi,
        R=np.ones((len(rows), 1)), doses=np.zeros((1, len(base))), offset=np.zeros(len(base))) for rows in folds]
    result = {}
    for family in ('NB', 'pooled'):
        normal = [profiles[family+'_'+label] for label in ('even', 'odd')]

        def score(parameter):
            value, error = 0., 0.
            for held in (0, 1):
                train = 1-held
                t = fraction(parameter, infos[train], mean_info)
                gain, _, rounding = helper.gain_gradient(validation[held], t*normal[train])
                value += gain
                error += rounding
            return value, error

        eps = np.finfo(float).eps
        coarse = minimize_scalar(lambda p: -score(p)[0], bounds=(0., 1.), method='bounded',
            options=dict(xatol=eps**.5, maxiter=500))
        fine = minimize_scalar(lambda p: -score(p)[0], bounds=(0., 1.), method='bounded',
            options=dict(xatol=eps**.75, maxiter=500))
        if not coarse.success or not fine.success:
            raise ArithmeticError('Information-scaled selection did not converge')
        candidates = [(0., *score(0.)), (1., *score(1.)), (float(fine.x), *score(fine.x))]
        parameter, value, error = max(candidates, key=lambda v: v[1])
        coarse_gain, coarse_error = score(coarse.x)
        assert value >= coarse_gain-error-coarse_error
        result[family] = dict(cv_parameter=parameter, t=fraction(parameter, full_info, mean_info),
            training_fractions=[fraction(parameter, v, mean_info) for v in infos],
            heldout_gain=value, numerical_budget=error, endpoint_zero_gain=candidates[0][1],
            endpoint_one_gain=candidates[1][1], coarse_parameter=float(coarse.x), fine_parameter=float(fine.x),
            normal_prediction_refinement=float(np.max(np.abs((fraction(coarse.x, full_info, mean_info)-
                fraction(fine.x, full_info, mean_info))*profiles[family+'_full']))),
            prior_precision=None if parameter == 0 else float(mean_info*(1-parameter)/parameter),
            zero_correction=parameter == 0, fold_information=infos, mean_fold_information=mean_info,
            full_information=full_info, full_information_additivity_budget=budget,
            scope='Trace-isotropic local approximation, not an exact NB posterior')
    return result
