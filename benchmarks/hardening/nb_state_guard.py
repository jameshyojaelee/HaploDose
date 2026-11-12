"""Research-only regional count comparison against continuous expression state."""
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logsumexp, roots_hermitenorm


def objective(x, O, E, r, order, mixture):
    """NB likelihood ratios integrated over one shared per-cell log-rate state.

    The omitted NB log-pmf at eta=0 is constant across both model families.
    Sigma is the latent log-rate SD, shared by both alternative components.
    """
    if x[-1] == 0:
        # At zero SD there is no integral. Evaluating its analytic limit also
        # preserves the exactly zero SD derivative instead of quadrature noise.
        means = x[:2] if mixture else x[:1]
        logs, scores = [], []
        baseline = np.logaddexp(np.log(r), np.log(E))
        for mean in means:
            shifted = np.log(E)+mean
            logs.append((O*mean-(O+r)*(np.logaddexp(np.log(r), shifted)-baseline)).sum(axis=1))
            scores.append((O-(O+r)*expit(shifted-np.log(r))).sum(axis=1))
        if mixture:
            joint = np.asarray(logs).T+np.array([-np.logaddexp(0.,x[2]),-np.logaddexp(0.,-x[2])])
            total = logsumexp(joint, axis=1)
            responsibility = np.exp(joint-total[:, None])
            gradient = np.r_[(responsibility*np.asarray(scores).T).sum(axis=0),
                (responsibility[:,1]-expit(x[2])).sum(), 0.]
        else:
            total = logs[0]
            gradient = np.array([np.sum(scores[0]), 0.])
        return -float(np.mean(total)), -gradient/len(O)
    nodes, weights = roots_hermitenorm(order)
    weights = np.log(weights)-.5*np.log(2*np.pi)
    means, sigma = (x[:2], x[3]) if mixture else (x[:1], x[1])
    logs, scores, sd_scores = [], [], []
    logE, logr = np.log(E), np.log(r)
    baseline = np.logaddexp(logr, logE)
    for mean in means:
        eta = mean+sigma*nodes
        shift = logE[:, :, None]+eta
        terms = O[:, :, None]*eta-(O+r)[:, :, None]*(np.logaddexp(logr[:, :, None], shift)-baseline[:, :, None])
        derivative = (O[:, :, None]-(O+r)[:, :, None]*expit(shift-logr[:, :, None])).sum(axis=1)
        joint = terms.sum(axis=1)+weights
        marginal = logsumexp(joint, axis=1)
        posterior = np.exp(joint-marginal[:, None])
        logs.append(marginal)
        scores.append((posterior*derivative).sum(axis=1))
        sd_scores.append((posterior*derivative*nodes).sum(axis=1))
    if mixture:
        logpi = np.array([-np.logaddexp(0., x[2]), -np.logaddexp(0., -x[2])])
        joint = np.array(logs).T+logpi
        total = logsumexp(joint, axis=1)
        R = np.exp(joint-total[:, None])
        gradient = np.r_[(R*np.array(scores).T).sum(axis=0),
                         (R[:, 1]-expit(x[2])).sum(),
                         (R*np.array(sd_scores).T).sum()]
    else:
        total = logs[0]
        gradient = np.array([np.sum(scores[0]), np.sum(sd_scores[0])])
    return -float(np.mean(total)), -gradient/len(O)


def variance_objective(x, O, E, r, order, mixture):
    """Optimize in variance, whose derivative need not vanish at its boundary."""
    x = np.asarray(x, float)
    if x[-1] < 0:
        raise ValueError('Nonnegative shared-state variance required')
    standard = x.copy()
    standard[-1] = np.sqrt(x[-1])
    value, gradient = objective(standard, O, E, r, order, mixture)
    if x[-1] > 0:
        gradient[-1] /= 2*standard[-1]
        return value, gradient
    # Gaussian heat equation: d_v E[f(mu+sqrt(v)Z)] at0 = f''(mu)/2.
    # f''/f = score_eta**2 + Hessian_eta. This removes the artificial
    # zero-SD stationary point without imposing a positive variance floor.
    means = x[:2] if mixture else x[:1]
    logs, derivatives = [], []
    for mean in means:
        shifted = np.log(E)+mean
        p = expit(shifted-np.log(r))
        score = (O-(O+r)*p).sum(axis=1)
        hessian = -((O+r)*p*(1-p)).sum(axis=1)
        derivatives.append(.5*(score**2+hessian))
        logs.append((O*mean-(O+r)*(np.logaddexp(np.log(r), shifted)-np.logaddexp(np.log(r), np.log(E)))).sum(axis=1))
    if mixture:
        joint = np.asarray(logs).T+np.array([-np.logaddexp(0., x[2]), -np.logaddexp(0., -x[2])])
        responsibility = np.exp(joint-logsumexp(joint, axis=1)[:, None])
        gradient[-1] = -float(np.mean((responsibility*np.asarray(derivatives).T).sum(axis=1)))
    else:
        gradient[-1] = -float(np.mean(derivatives[0]))
    return value, gradient


def fit_comparison(O, E, r, g, order):
    O, E, r = (np.asarray(v, float) for v in (O, E, r))
    r = np.broadcast_to(r, O.shape)
    if (O.ndim != 2 or O.shape != E.shape or not np.isfinite(O).all()
            or not np.isfinite(E).all() or not np.isfinite(r).all()
            or np.any(O < 0) or np.any(E <= 0) or np.any(r <= 0)):
        raise ValueError('Finite counts, supported means and positive NB sizes required')
    z = np.log((O.sum(axis=1)+.5)/(E.sum(axis=1)+.5))
    centre, spread = float(np.median(z)), float(np.std(z))
    starts_null = [[centre, 0.], [centre, spread**2]]
    options = dict(ftol=np.finfo(float).eps, gtol=np.sqrt(np.finfo(float).eps), maxiter=1000, maxls=50)
    nulls = [minimize(variance_objective, x, args=(O, E, r, order, False), jac=True,
                      bounds=[(None, None), (0., None)], method='L-BFGS-B', options=options) for x in starts_null]
    viable_null = [f for f in nulls if f.success and np.isfinite(f.fun)]
    if not viable_null:
        raise ArithmeticError('Continuous-state NB null optimization unresolved')
    null = min(viable_null, key=lambda f: f.fun)
    logit_bound = np.log(999.)  # Same proposal mixture-weight boundary as native code.
    means = [[float(z[g == k].mean()) for k in (0, 1)], list(np.quantile(z, [.2, .8]))]
    starts = [m+[0., s] for m in means for s in (0., float(null.x[1]))]
    starts.append([null.x[0], null.x[0], 0., null.x[1]])  # Nested null is feasible.
    alternatives = [minimize(variance_objective, x, args=(O, E, r, order, True), jac=True,
        bounds=[(None, None), (None, None), (-logit_bound, logit_bound), (0., None)],
        method='L-BFGS-B', options=options) for x in starts]
    viable = [f for f in alternatives if f.success and np.isfinite(f.fun)]
    if not viable:
        raise ArithmeticError('Discrete-plus-state NB optimization unresolved')
    alternative = min(viable, key=lambda f: f.fun)
    advantage = 2*len(O)*(null.fun-alternative.fun)-2*np.log(len(O))
    null_sd, alternative_sd = null.x.copy(), alternative.x.copy()
    null_sd[-1], alternative_sd[-1] = np.sqrt(null.x[-1]), np.sqrt(alternative.x[-1])
    return dict(order=order, BIC_advantage=float(advantage), null=null_sd.tolist(),
                alternative=alternative_sd.tolist(), null_objective=float(null.fun),
                alternative_objective=float(alternative.fun),
                optimized_state_parameter='variance',
                null_solutions=len(viable_null), alternative_solutions=len(viable))


def make_guard(records):
    def guard(self, cells, g, source, item):
        studies = []
        for order in (24, 48, 96, 192):
            result = fit_comparison(item['O'], item['E'], item['r'], g, order)
            studies.append(result)
            if len(studies) >= 2:
                discrepancy = abs(studies[-1]['BIC_advantage']-studies[-2]['BIC_advantage'])
                # Two successive quadrature orders must resolve the model choice.
                # This is a convergence diagnostic, not a formal integration bound.
                if abs(result['BIC_advantage']) > 2*discrepancy:
                    break
        resolved = abs(result['BIC_advantage']) > 2*discrepancy
        passed = bool(resolved and result['BIC_advantage'] > 0)
        records.append(dict(source=source, proposal_cells=np.asarray(cells).tolist(),
            specimen_rows=self.rows.rows.tolist(), quadrature=studies,
            quadrature_choice_resolved=bool(resolved), passed=passed,
            interpretation='Conditional working-model BIC; not calibrated or genomic evidence'))
        return passed
    return guard
