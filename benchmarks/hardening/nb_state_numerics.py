"""Independent numerical references for the regional count-state safeguard."""
import argparse
import importlib.util
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit
from scipy.stats import nbinom, norm

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    path = Path(__file__).with_name('nb_state_guard.py')
    spec = importlib.util.spec_from_file_location('checked_nb_state', path)
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    run = start_run(args.out, dict(script=Path(__file__), model=path, design=args.task/'DESIGN.md'), {},
                    dict(scope='Numerical likelihood/gradient checks, not scientific calibration'), mode='nb_state_numerics')
    try:
        O = np.array([[2., 4., 1.], [7., 3., 5.], [1., 0., 2.], [4., 2., 3.]])
        E = np.array([[3., 2., 1.], [5., 4., 3.], [2., 1., 2.], [4., 3., 2.]])
        r = np.broadcast_to([3., 7., 11.], O.shape)
        base = nbinom.logpmf(O, r, r/(r+E)).sum(axis=1)
        rows = []
        for mixture, x in [(False, np.array([.1, .3])), (True, np.array([-.2, .4, -.3, .3])),
                           (False, np.array([.1, 0.])), (True, np.array([.1, .1, -.7, .3]))]:
            means, sigma = (x[:2], x[3]) if mixture else (x[:1], x[1])
            weights = [1-expit(x[2]), expit(x[2])] if mixture else [1.]
            logvalues, allowances = [], []
            for y, e, size in zip(O, E, r):
                probability, error = 0., 0.
                for mean, weight in zip(means, weights):
                    def integrand(z):
                        mu = e*np.exp(mean+sigma*z)
                        return float(np.exp(nbinom.logpmf(y, size, size/(size+mu)).sum())*norm.pdf(z))
                    # The omitted normal tail is bounded by 2*sf(12), since
                    # a product of count probabilities cannot exceed one.
                    integral, estimate = quad(integrand, -12., 12., epsabs=1e-13, epsrel=1e-13)
                    probability += weight*integral
                    error += weight*(estimate+2*norm.sf(12.))
                assert probability > error
                logvalues.append(np.log(probability))
                allowances.append(error/(probability-error))
            expected = -float(np.mean(np.array(logvalues)-base))
            values = [model.objective(x, O, E, r, order, mixture)[0] for order in (48, 96, 192)]
            rounding = 256*O.size*np.finfo(float).eps*(1+np.max(np.abs(base))+abs(expected))
            allowance = float(np.mean(allowances))+rounding
            assert abs(values[-1]-expected) <= allowance, (values, expected, allowance)
            value, gradient = model.objective(x, O, E, r, 96, mixture)
            # Independent central differences, with error assessed from two
            # step sizes and the cancellation scale of these finite objectives.
            differences = []
            for j in range(len(x)):
                step = np.cbrt(np.finfo(float).eps)*(1+abs(x[j]))
                estimates = []
                for h in (step, step/2):
                    plus, minus = x.copy(), x.copy()
                    plus[j] += h; minus[j] -= h
                    estimates.append((model.objective(plus,O,E,r,96,mixture)[0]-
                                      model.objective(minus,O,E,r,96,mixture)[0])/(2*h))
                bound = 2*abs(estimates[1]-estimates[0])+rounding/step
                assert abs(gradient[j]-estimates[-1]) <= bound, (gradient[j], estimates, bound)
                differences.append(abs(float(gradient[j]-estimates[-1])))
            rows.append(dict(mixture=mixture, parameters=x.tolist(), objective=value,
                independent_scipy_quad=expected, discrepancy=abs(value-expected),
                quad_plus_rounding_allowance=allowance, gradient_discrepancies=differences,
                quadrature_convergence=values))
        null = np.array([.2, .4])
        alt = np.array([.2, .2, .7, .4])
        nested = model.objective(null,O,E,r,96,False)[0]
        mixture = model.objective(alt,O,E,r,96,True)[0]
        assert abs(nested-mixture) <= rounding
        write_json(dict(checks=rows, identical_components_reproduce_null=True,
                        all_checks_passed=True), args.out/'report.json')
        finish_run(args.out,run,backend='CPU independent SciPy NB PMF/adaptive quadrature/finite differences')
    except BaseException as error:
        finish_run(args.out,run,error=error)
        raise


if __name__ == '__main__':
    main()
