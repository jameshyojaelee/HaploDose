"""Independent Fisher/quadratic checks and exact904 default replay before905 data."""
import argparse
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.optimize import minimize_scalar
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    info = load('info', root/'normal_offset_information.py')
    selection = load('selection', root/'normal_offset_shrink.py')
    helper = load('objective', root/'offset_nb_objective.py')
    estimator = load('estimator', root/'reference_offset_panel.py')
    prior = a.task.parent/'normal_offset_shrink_20261109'
    rec = start_run(a.out, dict(script=Path(__file__), info=root/'normal_offset_information.py',
        selector=root/'normal_offset_shrink.py', generator=root/'normal_offset_prepare.py',
        design=a.task/'DESIGN.md', prior_observations=prior/'panel_v1/observations.npz',
        prior_reference=prior/'fit_v1/references.json'), {},
        dict(scope='Numerical gate and retrospective default regression, biological n=0'), mode='normal_information_numerics')
    try:
        eps = np.finfo(float).eps
        base, phi = np.array([.1, .2, .3, .4]), np.array([.1, .2, .05, .3])
        L = np.array([1200., 2400., 600., 4800.])
        nuisance = np.log(np.array([[1., 2., 1., .5], [2., 1., .5, 1.], [1., 1., 1., 1.], [.5, 2., 2., 1.]]))
        traces = []
        for h in (eps**.75, eps**.75/2):
            total = np.zeros((4, 4))
            for depth, row in zip(L, nuisance):
                logits = np.log(base)+row

                def logmean(z):
                    rates = np.exp(z)
                    return np.log(depth*rates/sum(rates))

                derivative = np.column_stack([np.imag(logmean(logits.astype(complex)+
                    1j*h*np.eye(4)[j]))/h for j in range(4)])
                mean = np.exp(logmean(logits))
                # E(score^2)=Var(Y)/(1+phi*m)^2 derives native NB weights.
                weight = (mean+phi*mean**2)/(1+phi*mean)**2
                total += derivative.T@(weight[:, None]*derivative)
            traces.append(float(np.trace(total)/3))
        fisher_value = info.information(L, base, nuisance, phi)
        budget = 128*L.size*base.size*eps*max(1., fisher_value)
        assert max(abs(v-fisher_value) for v in traces) <= budget
        assert abs(sum(info.information(L[p::2], base, nuisance[p::2], phi) for p in (0, 1))-fisher_value) <= budget
        assert abs(info.information(np.tile(L, 2), base, np.tile(nuisance, (2, 1)), phi)-2*fisher_value) <= 2*budget
        assert abs(info.information(L, base*np.exp(7), nuisance, phi)-fisher_value) <= budget
        assert abs(info.information(L, base[::-1], nuisance[:, ::-1], phi[::-1])-fisher_value) <= budget
        quadratic = []
        for lam in (0, 1, 3, 9):
            parameter = Fraction(4, 4+lam)
            for I in (2, 6, 8):
                actual = info.fraction(float(parameter), I, 4)
                expected = Fraction(I, I+lam)
                assert abs(actual-float(expected)) <= 128*eps
                y = Fraction(7, 10)
                solution = expected*y
                objective = lambda x: .5*I*(x-float(y))**2+.5*lam*x*x
                result = minimize_scalar(objective, method='brent', options=dict(xtol=eps**.75, maxiter=500))
                assert result.success
                rounding = 128*eps*max(1., objective(float(solution)), objective(result.x))
                assert abs(objective(result.x)-objective(float(solution))) <= rounding
                position_budget = np.sqrt(2*rounding/(I+lam))
                assert abs(result.x-float(solution)) <= position_budget
                quadratic.append(dict(prior_precision=lam, information=I, exact_fraction=str(expected),
                    computed_fraction=actual, scipy_solution=float(result.x), exact_solution=str(solution),
                    position_budget=float(position_budget)))
        for parameter in (0., .25, .5, .75, 1.):
            assert abs(info.fraction(parameter, 8, 4)-2*parameter/(1+parameter)) <= 128*eps
        H, lam, y = np.diag([1., 4., 9.]), 2., np.array([-.4, .2, .7])
        exact = np.linalg.solve(H+lam*np.eye(3), H@y)
        average = np.trace(H)/3
        scalar = average/(average+lam)*y
        objective = lambda x: .5*(x-y)@H@(x-y)+.5*lam*(x@x)
        assert objective(scalar) > objective(exact)

        # Recreate the unmodified default seed and independently bind every array.
        subprocess.run([sys.executable, '-I', str(root/'normal_offset_prepare.py'), '--out',
            str(a.out/'default904'), '--task', str(a.task)], check=True)
        array_checks = 0
        for name in ('observations.npz', 'targets.npz', 'observations_query_changed.npz', 'targets_changed.npz'):
            old, new = np.load(prior/'panel_v1'/name), np.load(a.out/'default904'/name)
            assert set(old.files) == set(new.files)
            for key in old.files:
                assert old[key].dtype == new[key].dtype and old[key].shape == new[key].shape
                assert old[key].tobytes() == new[key].tobytes()
                array_checks += 1
        assert (prior/'panel_v1/metadata.json').read_bytes() == (a.out/'default904/metadata.json').read_bytes()
        data, old_arrays = np.load(prior/'panel_v1/observations.npz'), np.load(prior/'fit_v1/reference_arrays.npz')
        old_references = json.loads((prior/'fit_v1/references.json').read_text())
        replayed = 0
        for key, expected in old_references.items():
            fields = {name: data[key+'__'+name] for name in ('N', 'Ln', 'base', 'nn', 'phi')}
            observed, arrays = selection.normal_selection(fields['N'], fields['Ln'], fields['base'], fields['nn'],
                fields['phi'], helper, estimator)
            assert observed == expected
            for name, value in arrays.items():
                assert value.tobytes() == old_arrays[key+'__'+name].tobytes()
                replayed += 1
        (a.out/'report.json').write_text(json.dumps(dict(all_checks_passed=True, Fisher_information=fisher_value,
            independent_complex_step_traces=traces, Fisher_budget=budget, quadratic=quadratic,
            anisotropic_counterexample=dict(exact=exact.tolist(), scalar=scalar.tolist(),
                exact_objective=float(objective(exact)), scalar_objective=float(objective(scalar)),
                scope='Scalar information is approximate outside the isotropic quadratic limit'),
            default904_arrays_bit_identical=array_checks, default904_normal_arrays_bit_identical=replayed,
            default904_reference_records_identical=len(old_references), biological_specimens=0,
            production_changed=False, scientific_criteria_changed=False), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU complex-step Fisher, exact quadratic/SciPy and default904 bit replay')


if __name__ == '__main__':
    main()
