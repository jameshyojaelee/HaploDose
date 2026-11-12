"""Independent scalar NB reference and Decimal80 audit of conditional offsets."""
import argparse
from decimal import Decimal, localcontext
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import softmax
from scipy.stats import nbinom
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def full_q(data, offset, doses=None):
    doses = data['doses'] if doses is None else doses
    out = 0.
    for k, dose in enumerate(doses):
        mean = data['L'][:, None]*softmax(np.log(data['base'])+offset+dose+data['nuisance'], axis=1)
        size = 1/data['phi']
        value = nbinom.logpmf(data['Y'], size, size/(size+mean))
        out += float(np.sum(data['R'][:, k, None]*value))
    return out


def decimal_gain(data, offset, doses=None):
    """Independent scalar probabilities and NB ratios, avoiding softmax/logaddexp."""
    doses = data['doses'] if doses is None else doses
    D = lambda value: Decimal(str(float(value)))
    with localcontext() as context:
        context.prec = 80
        total = Decimal(0)
        for i, counts in enumerate(data['Y']):
            for k in range(data['R'].shape[1]):
                if data['R'][i, k] == 0:
                    continue
                rates = []
                for o, ds in ((data['offset'], data['doses']), (offset, doses)):
                    raw = [D(data['base'][j])*(D(o[j])+D(ds[k, j])+D(data['nuisance'][i, j])).exp()
                           for j in range(len(counts))]
                    rates.append([D(data['L'][i])*v/sum(raw) for v in raw])
                for j, y in enumerate(counts):
                    m0, m = rates[0][j], rates[1][j]
                    r = 1/D(data['phi'][j])
                    total += D(data['R'][i, k])*(D(y)*(m/m0).ln()-(D(y)+r)*((r+m)/(r+m0)).ln())
        return float(total)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    helper = load('nb_offset', root/'offset_nb_objective.py')
    fixture = load('counterexample', root/'offset_expectation.py')
    conditional = load('conditional', root/'offset_conditional.py')
    rec = start_run(a.out, dict(script=Path(__file__), helper=root/'offset_nb_objective.py',
        fixture=root/'offset_expectation.py', design=a.task/'DESIGN.md'), {},
        dict(scope='Fixed conditional two-bin numerical reference, biological n=0'), mode='offset_nb_numerics')
    try:
        cases = []
        eps = np.finfo(float).eps
        for condition in ('none', 'depth', 'programme', 'both'):
            for membership in ('hard', 'dropped', 'soft'):
                for signal in (False, True):
                    model, R, doses, keep, _, _ = fixture.fixture(condition, membership, signal, physical_depth=True)
                    model.phi = np.array([.2, .35])
                    data = helper.arrays(model, R[:, keep], doses)
                    fitted, diagnostics = helper.optimize(data)
                    conditional.update_offset(model, R, doses, keep)
                    gain, _, budget = helper.gain_gradient(data, fitted)
                    reference = decimal_gain(data, fitted)
                    assert abs(gain-reference) <= budget
                    # Search bounds locate a reference bracket, never constrain the fit.
                    grid = np.linspace(-16., 16., 257)
                    values = [-full_q(data, np.array([x, 0.])) for x in grid]
                    index = int(np.argmin(values))
                    assert 0 < index < len(grid)-1
                    scalar = minimize_scalar(lambda z: -full_q(data, np.array([z, 0.])),
                        bracket=grid[index-1:index+2], method='brent', options=dict(xtol=eps**.75, maxiter=500))
                    assert scalar.success
                    scalar_offset = np.array([scalar.x/2, -scalar.x/2])
                    scalar_gain = decimal_gain(data, scalar_offset)
                    assert abs(gain-scalar_gain) <= budget
                    trial = np.array([.13, -.13])
                    analytic = helper.gain_gradient(data, trial)[1][0]
                    h = eps**(1/3)
                    fd = []
                    for step in (h, h/2, h/4):
                        plus, minus = trial.copy(), trial.copy()
                        plus[0] += step
                        minus[0] -= step
                        fd.append((full_q(data, plus)-full_q(data, minus))/(2*step))
                    # Refinement discrepancies plus full-logpmf subtraction roundoff.
                    fd_budget = 8*max(abs(fd[1]-fd[0]), abs(fd[2]-fd[1]))+128*eps*max(1., abs(full_q(data, trial)))/(h/4)
                    assert abs(analytic-fd[-1]) <= fd_budget
                    gauge_gain = helper.gain_gradient(data, fitted+7.)[0]
                    assert abs(gauge_gain-gain) <= budget
                    conditional_gain = decimal_gain(data, model.offset)
                    assert gain >= conditional_gain-budget
                    cases.append(dict(condition=condition, membership=membership, signal=signal,
                        fitted_offset=fitted.tolist(), scalar_offset=scalar_offset.tolist(),
                        gain=gain, Decimal_gain=reference, scalar_Decimal_gain=scalar_gain,
                        conditional_gain=conditional_gain, budget=budget, analytic_gradient=float(analytic),
                        finite_differences=fd, finite_difference_budget=fd_budget, optimizer=diagnostics))
        (a.out/'report.json').write_text(json.dumps(dict(cases=cases, all_checks_passed=True,
            Decimal_precision=80, production_changed=False, biological_specimens=0,
            scientific_criteria_changed=False), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU SciPy NB logpmf/Brent, finite differences, scalar Decimal80')


if __name__ == '__main__':
    main()
