"""Exact derivative and missed-null counterexamples for variance optimization."""
import argparse
from fractions import Fraction
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    path = Path(__file__).with_name('nb_state_guard.py')
    spec = importlib.util.spec_from_file_location('checked_variance_state', path)
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    rec = start_run(args.out, dict(script=Path(__file__), model=path,
        design=args.task/'CORRECTION_DESIGN.md', counterexamples=args.task/'boundary_v2/report.json'), {},
        dict(scope='Optimization correctness; no sampling or biological claim'), mode='nb_state_boundary_regression')
    try:
        O = np.array([[0.], [0.], [8.], [8.]])
        E, r = np.full_like(O, 4.), np.full_like(O, 1000.)
        p = Fraction(4, 1004)
        exact = -sum(((Fraction(int(y))-Fraction(int(y)+1000)*p)**2
            -Fraction(int(y)+1000)*p*(1-p))/2 for y in O.ravel())/4
        value, gradient = model.variance_objective(np.array([0., 0.]), O, E, r, 96, False)
        operations = 256*O.size+256
        bound = operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)*max(1., abs(float(exact)))
        assert abs(gradient[-1]-float(exact)) <= bound and gradient[-1] < 0
        assert model.objective(np.array([0., 0.]), O, E, r, 96, False)[1][-1] == 0
        grid = []
        for v in (1e-3, 5e-4, 1e-4, 5e-5):
            forward = (model.variance_objective(np.array([0., v]), O, E, r, 96, False)[0]-value)/v
            grid.append(dict(variance=v, forward_derivative=forward))
        assert abs(grid[-1]['forward_derivative']-float(exact)) < abs(grid[0]['forward_derivative']-float(exact))
        fit = minimize(model.variance_objective, [0., 0.], args=(O,E,r,96,False), jac=True,
            bounds=[(None,None),(0.,None)], method='L-BFGS-B',
            options=dict(ftol=np.finfo(float).eps, gtol=np.sqrt(np.finfo(float).eps), maxiter=1000, maxls=50))
        assert fit.success and fit.x[-1] > 0 and fit.fun < value-bound
        # Both null and mixture boundary derivatives must match an independent
        # change in the original SD parameter at a strictly positive variance.
        interior = []
        for mixture, x in [(False,np.array([.1,.04])),(True,np.array([-.2,.3,-.7,.04]))]:
            f,g = model.variance_objective(x,O,E,r,96,mixture)
            errors = []
            for j in range(len(x)):
                h = np.cbrt(np.finfo(float).eps)*(1+abs(x[j]))
                values = []
                for step in (h,h/2):
                    plus,minus=x.copy(),x.copy();plus[j]+=step;minus[j]-=step
                    values.append((model.variance_objective(plus,O,E,r,96,mixture)[0]-
                        model.variance_objective(minus,O,E,r,96,mixture)[0])/(2*step))
                allowance=2*abs(values[1]-values[0])+bound/h
                assert abs(g[j]-values[-1]) <= allowance
                errors.append(float(abs(g[j]-values[-1])))
            interior.append(dict(mixture=mixture, gradient_errors=errors))
        write_json(dict(exact_variance_boundary_derivative=str(exact), computed=float(gradient[-1]),
            boundary_error=float(abs(gradient[-1]-float(exact))), rounding_budget=bound,
            sd_zero_stationary_point_reproduced=True, variance_start_zero_escapes=True,
            fitted_variance=float(fit.x[-1]), objective_decrease=float(value-fit.fun),
            forward_study=grid, interior=interior, no_positive_floor=True), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU exact Fraction variance derivative and optimization counterexample')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':main()
