"""Research-only conditional NB objective; no DNA or joint-EM optimality claim."""
from contextlib import contextmanager

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from haplodose.clones import SpecimenModel


def arrays(model, R, doses, offset=None):
    Y = np.concatenate([np.asarray(y, float) for _, _, y in model.rows.chunks()])
    return dict(Y=Y, L=np.asarray(model.L, float).copy(), R=np.asarray(R, float).copy(),
        doses=np.asarray(doses, float).copy(), base=np.asarray(model.base, float).copy(),
        nuisance=(model.x[:, None]*model.gamma[None]+model.U@model.W).copy(),
        phi=np.asarray(model.phi, float).copy(),
        offset=np.asarray(model.offset if offset is None else offset, float).copy())


def gain_gradient(data, offset, doses=None):
    """NB Q gain against stored offset/doses; analytic gradient w.r.t. shared logits.

    Holding latent weights fixed, Q=sum_ik R_ik sum_j log NB(Y_ij;m_ikj,phi_j),
    m=L softmax(log base+offset+dose_k+nuisance_i). Terms constant in m cancel.
    Accumulates n x J temporaries per clone, never an n x K x J tensor.
    """
    Y, L, R = data['Y'], data['L'], data['R']
    doses = data['doses'] if doses is None else np.asarray(doses, float)
    value, magnitude, gradient = 0., 0., np.zeros(Y.shape[1])
    logbase, invphi = np.log(data['base']), 1/data['phi']
    for k in range(R.shape[1]):
        eta0 = logbase+data['offset']+data['doses'][k]+data['nuisance']
        lp0 = eta0-logsumexp(eta0, axis=1, keepdims=True)
        p0 = np.exp(lp0)
        delta = (offset-data['offset'])+(doses[k]-data['doses'][k])
        delta -= np.median(delta)  # an exact composition gauge; avoid common large logits
        d = delta-np.log1p(np.sum(p0*np.expm1(delta), axis=1, keepdims=True))
        probabilities = p0*np.exp(d)
        m0 = L[:, None]*p0
        m = L[:, None]*probabilities
        a = Y*d
        b = (Y+invphi)*np.log1p(m0/(invphi+m0)*np.expm1(d))
        w = R[:, k, None]
        value += float(np.sum(w*(a-b)))
        magnitude += float(np.sum(w*(np.abs(a)+np.abs(b))))
        score = (Y-m)/(1+data['phi']*m)
        gradient += np.sum(w*(score-probabilities*score.sum(axis=1, keepdims=True)), axis=0)
    ops = 128*Y.size*max(1, R.shape[1])
    eps = np.finfo(float).eps
    budget = ops*eps/(1-ops*eps)*max(1., magnitude)
    return value, gradient, budget


def optimize(data):
    """Gauge-fixed L-BFGS-B, with explicit convergence and objective diagnostics."""
    start = data['offset']-data['offset'][-1]
    scale = max(1., float(np.sum(data['R'].sum(axis=1)[:, None]*data['Y'])))

    def objective(z):
        value, gradient, _ = gain_gradient(data, np.r_[z, 0.])
        return -value/scale, -gradient[:-1]/scale

    eps = np.finfo(float).eps
    result = minimize(objective, start[:-1], method='L-BFGS-B', jac=True,
        options=dict(gtol=eps**.75, ftol=16*eps, maxiter=500, maxls=40))
    fitted = np.r_[result.x, 0.]
    fitted -= np.median(fitted)
    gain, gradient, budget = gain_gradient(data, fitted)
    if not np.isfinite(fitted).all() or gain < -budget:
        raise FloatingPointError('Conditional optimizer lost its anchored NB objective')
    return fitted, dict(success=bool(result.success), message=str(result.message),
        iterations=int(result.nit), objective_gain=gain, budget=budget,
        max_abs_gradient=float(np.max(np.abs(gradient))), scale=scale,
        normalized_gradient=float(np.max(np.abs(gradient)))/scale)


@contextmanager
def objective_trace(arm, trace, worst, *, trace_applied=False):
    """Trace profile/offset updates; only conditional/NB arms alter offset updates."""
    from pathlib import Path
    import importlib.util
    spec = importlib.util.spec_from_file_location('conditional_offset', Path(__file__).with_name('offset_conditional.py'))
    conditional = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(conditional)
    names = ('em', 'pooled', 'loglik', 'shared_profiles', '_update_offset')
    originals = {name: getattr(SpecimenModel, name) for name in names}
    states, active = {}, set()

    def state(model):
        return states.setdefault(id(model), dict(model=len(states), instance=model))

    def em(model, *args, **kwargs):
        active.add(id(model))
        try:
            return originals['em'](model, *args, **kwargs)
        finally:
            active.discard(id(model))

    def pooled(model, R, *args, **kwargs):
        if id(model) in active:
            state(model)['R'] = np.asarray(R, float).copy()
        return originals['pooled'](model, R, *args, **kwargs)

    def loglik(model, doses, *args, **kwargs):
        if id(model) in active:
            state(model)['doses'] = np.asarray(doses, float).copy()
        return originals['loglik'](model, doses, *args, **kwargs)

    def record(model, kind, data, offset, doses=None, diagnostics=None):
        value, gradient, budget = gain_gradient(data, offset, doses)
        row = dict(index=len(trace), model=state(model)['model'], kind=kind, arm=arm,
            cells=model.rows.n, bins=model.J, components=data['R'].shape[1], gain=value,
            budget=budget, decreased=bool(value < -budget), max_abs_gradient=float(np.max(np.abs(gradient))))
        if diagnostics is not None:
            row['optimizer'] = diagnostics
        trace.append(row)
        key = (row['model'], kind)
        if key not in worst or value < worst[key]['row']['gain']:
            worst[key] = dict(row=row, data=data, candidate_offset=np.asarray(offset).copy(),
                candidate_doses=np.asarray(data['doses'] if doses is None else doses).copy())

    def profiles(model, *args, **kwargs):
        info = state(model)
        data = None
        if id(model) in active and 'R' in info and 'doses' in info and info['R'].shape[1] == len(info['doses']):
            data = arrays(model, info['R'], info['doses'])
        result = originals['shared_profiles'](model, *args, **kwargs)
        if data is not None:
            record(model, 'profile', data, model.offset, result['dose'])
            if trace_applied:
                info['profile_data'] = data
        return result

    def offset(model, R, doses, keep):
        if trace_applied:
            previous = state(model).pop('profile_data')
            assert np.array_equal(previous['R'], R)
            assert keep[0] == 0 and len(set(keep)) == len(keep)
            assert np.array_equal(np.asarray(doses)[0], np.zeros(model.J))
            assert np.array_equal(previous['offset'], model.offset)
            retained = dict(previous, R=previous['R'][:, keep], doses=previous['doses'][keep])
            record(model, 'profile_applied', retained, model.offset, doses)
            trace[-1].update(keep=list(keep), old_components=R.shape[1],
                responsibility_dtype=np.asarray(R).dtype.str,
                actual_normal_dose_zero=True, retained_mass_not_renormalized=True)
            item = worst[(trace[-1]['model'], 'profile_applied')]
            if item['row'] is trace[-1]:
                item['extra'] = dict(pre_keep_R=previous['R'], pre_keep_doses=previous['doses'],
                    keep=np.asarray(keep))
        data = arrays(model, np.asarray(R)[:, keep], doses)
        if arm == 'legacy':
            originals['_update_offset'](model, R, doses, keep)
            diagnostics = None
        elif arm == 'conditional':
            conditional.update_offset(model, R, doses, keep)
            diagnostics = None
        elif arm == 'NB':
            model.offset, diagnostics = optimize(data)
        else:
            raise ValueError(arm)
        record(model, 'offset', data, model.offset, diagnostics=diagnostics)

    for name, method in zip(names, (em, pooled, loglik, profiles, offset)):
        setattr(SpecimenModel, name, method)
    try:
        yield
    finally:
        for name, method in originals.items():
            setattr(SpecimenModel, name, method)
