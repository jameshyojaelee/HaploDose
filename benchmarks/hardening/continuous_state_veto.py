"""Research guard: discrete regional RNA populations must beat continuous state."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logsumexp

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def load(name):
    spec = importlib.util.spec_from_file_location('continuous_guard_'+name, Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def compare(z, measurement_variance, g):
    z, measurement_variance, g = map(np.asarray, (z, measurement_variance, g))
    assert z.shape == measurement_variance.shape == g.shape
    assert np.isfinite(z).all() and np.isfinite(measurement_variance).all() and np.all(measurement_variance > 0)
    scale = float(np.std(z))
    if scale == 0:
        return dict(passed=False, reason='no_RNA_variation')
    centre = float(np.mean(z))
    z, v = (z-centre)/scale, measurement_variance/scale**2
    n = len(z)
    # One Gaussian latent regional state, convolved with known working errors.
    def null(x):
        V = v+x[1]
        residual = z-x[0]
        ll = -.5*(np.log(2*np.pi*V)+residual**2/V)
        grad = np.array([(residual/V).sum(), .5*(residual**2/V**2-1/V).sum()])
        return -ll.mean(), -grad/n
    eps = np.finfo(float).eps
    options = dict(ftol=eps, gtol=np.sqrt(eps), maxiter=1000, maxls=50)
    fitted_null = minimize(null, [0., max(0., 1-float(v.mean()))], jac=True,
                           method='L-BFGS-B', bounds=[(None, None), (0., None)], options=options)
    # Two discrete locations with the SAME continuous state variance in both.
    # Shared variance prevents a small component collapsing onto one cell.
    def alternative(x):
        V = v+x[3]
        residual = z[:, None]-x[None, :2]
        pi = expit(x[2])
        ll = -.5*(np.log(2*np.pi*V[:, None])+residual**2/V[:, None])
        ll += [-np.logaddexp(0., x[2]), -np.logaddexp(0., -x[2])]
        total = logsumexp(ll, axis=1)
        R = np.exp(ll-total[:, None])
        grad = np.r_[(R*residual/V[:, None]).sum(axis=0), (R[:, 1]-pi).sum(),
                     .5*(R*(residual**2/V[:, None]**2-1/V[:, None])).sum()]
        return -total.mean(), -grad/n
    starts = [[float(z[g == k].mean()) for k in (0, 1)], list(np.quantile(z, [.2, .8]))]
    fits = [minimize(alternative, means+[0., max(0., float(fitted_null.x[1])/2)], jac=True,
                    method='L-BFGS-B', bounds=[(None, None), (None, None),
                    (-np.log(999.), np.log(999.)), (0., None)], options=options) for means in starts]
    viable = [f for f in fits if f.success and np.isfinite(f.fun) and np.isfinite(f.x).all()]
    if not fitted_null.success or not viable:
        return dict(passed=False, reason='optimization_unresolved', null_message=str(fitted_null.message),
                    alternative_messages=[str(f.message) for f in fits])
    best = min(viable, key=lambda f: f.fun)
    gain = float(n*(fitted_null.fun-best.fun))
    bic_difference = 2*gain-2*np.log(n)
    from scipy.stats import norm
    null_ll = float(norm.logpdf(z, fitted_null.x[0], np.sqrt(v+fitted_null.x[1])).sum())
    pi = expit(best.x[2])
    alternative_ll = float(np.logaddexp(
        norm.logpdf(z, best.x[0], np.sqrt(v+best.x[3]))+np.log1p(-pi),
        norm.logpdf(z, best.x[1], np.sqrt(v+best.x[3]))+np.log(pi)).sum())
    operations = 128*n+64
    rounding = operations*eps/(1-operations*eps)*(1+abs(null_ll)+abs(alternative_ll))
    if abs(null_ll+n*fitted_null.fun) > rounding or abs(alternative_ll+n*best.fun) > rounding:
        raise ArithmeticError('Independent Normal likelihood disagrees with optimized objective')
    if abs(bic_difference) <= 4*rounding:
        return dict(passed=False, reason='numerically_unresolved_BIC', rounding_budget=rounding)
    return dict(passed=bool(bic_difference > 0), reason='working_BIC', cells=n,
        log_likelihood_gain=gain, discrete_BIC_advantage=float(bic_difference),
        null_parameters=fitted_null.x.tolist(), alternative_parameters=best.x.tolist(),
        centering=centre, scaling=scale, local_optimization_only=True,
        independent_likelihood_rounding_budget=float(rounding))


def make_guard(records):
    def guard(self, cells, g, source, item):
        O, E, r = item['O'], item['E'], item['r']
        assert O.shape == E.shape == r.shape and len(O) == len(cells)
        expected = E.sum(axis=1)
        observed = O.sum(axis=1)
        z = np.log((observed+.5)/(expected+.5))
        variance = (E+E*E/r).sum(axis=1)/(expected+.5)**2
        result = compare(z, variance, g)
        records.append(dict(source=source, specimen_rows=self.rows.rows.tolist(),
                            proposal_cells=np.asarray(cells).tolist(), **result))
        return result['passed']
    return guard


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
                  stage=Path(__file__).with_name('membership_weighted_stages.py'))
    rec = start_run(args.out, inputs, dict(BIC_extra_parameters=2),
        dict(scope='Approximate RNA-state model comparison, no biological or calibration guarantee'),
        mode='continuous_state_veto')
    records = []
    import sys
    argv = sys.argv
    try:
        noise = load('reference_gene_noise')
        sys.argv = [str(noise.__file__), '--out', str(args.out/'candidate'), '--task',
            str(args.task.parent/'depth_stratified_gene_noise_20261109'), '--stage', 'both', '--depth-strata', '5',
            '--continuous-state-guard']
        noise.main()
        records = json.loads((args.out/'candidate/panel/guard_records.json').read_text())
        report = json.loads((args.out/'candidate/report.json').read_text())
        write_json(dict(recovered=report['recovered'], cases=report['cases'],
                        state_guard_calls=len(records), production_changed=False), args.out/'report.json')
        write_json(records, args.out/'guard_records.json')
        finish_run(args.out, rec, backend='CPU complete recovery with continuous-state alternative guard')
    except BaseException as error:
        write_json(records, args.out/'guard_records.json')
        finish_run(args.out, rec, error=error)
        raise
    finally:
        sys.argv = argv


if __name__ == '__main__':
    main()
