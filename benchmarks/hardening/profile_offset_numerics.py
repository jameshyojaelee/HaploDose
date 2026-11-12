"""Independent normalized means for offset-conditioned native profile baselines."""
import argparse
from decimal import Decimal, localcontext
import importlib.util
import json
from pathlib import Path

import numpy as np
from haplodose.clones import CloneParams, Rows, SpecimenModel
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def independent_expected(base, offset, original_gamma, Lk, XLk):
    """Scalar Decimal80 depth-moment approximation, not exact depth integration."""
    D = lambda x: Decimal.from_float(float(x))
    with localcontext() as context:
        context.prec = 80
        probabilities = [D(b)*D(o).exp() for b, o in zip(base, offset)]
        total = sum(probabilities)
        probabilities = [p/total for p in probabilities]
        gamma = [D(g) for g in original_gamma]
        if np.any(offset != offset[0]):
            origin = sum(p*g for p, g in zip(probabilities, gamma))
            gamma = [g-origin for g in gamma]
        curvature = sum(p*g*g for p, g in zip(probabilities, gamma))
        result = []
        for depth, moments in zip(np.atleast_1d(Lk), np.atleast_2d(XLk)):
            mean = D(moments[0])/max(D(depth), Decimal(1))
            variance = max(Decimal(0), D(moments[1])/max(D(depth), Decimal(1))-mean*mean)
            raw = [p*((g*mean+g*g*variance/2)/(1+curvature*variance)).exp()
                   for p, g in zip(probabilities, gamma)]
            result.append([float(D(depth)*v/sum(raw)) for v in raw])
    return np.array(result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    hook_path = Path(__file__).with_name('profile_offset.py')
    hook = load('profile_hook', hook_path)
    rec = start_run(a.out, dict(script=Path(__file__), hook=hook_path, design=a.task/'DESIGN.md'), {},
        dict(scope='Known conditional RNA mean and depth-moment arithmetic, biological n=0'), mode='profile_offset_numerics')
    try:
        Y = np.tile(np.r_[np.full(20, 100), np.full(20, 200)], (4, 1))
        reset = np.zeros(40, bool)
        reset[0] = True
        model = SpecimenModel(Rows(Y, np.arange(4), 2), Y.sum(axis=1), np.full((40, 1), 1/40), reset,
            np.full(40, .2), CloneParams(n_programmes=0, spatial_neighbors=0))
        Lk = np.array([24000.])
        XLk = np.zeros((1, 2))
        counts = Y.sum(axis=0)[None]
        eps = np.finfo(float).eps
        model.offset = np.r_[np.zeros(20), np.full(20, np.log(2))]
        legacy = model.profiles(counts, Lk, XLk)
        trace, contexts = [], {}
        original_base, original_gamma = model.base.copy(), model.gamma.copy()
        with hook.conditioned_profiles(trace, contexts):
            corrected = model.profiles(counts, Lk, XLk)
        budget = 128*40*eps*float(Lk[0])
        assert np.max(np.abs(corrected['E']-counts)) <= budget
        assert np.max(np.abs(legacy['E']-counts)) > budget
        assert np.array_equal(model.base, original_base) and np.array_equal(model.gamma, original_gamma)
        cases = []
        for baseline in ('flat', 'heterogeneous'):
            base = np.ones(40) if baseline == 'flat' else np.arange(1, 41, dtype=float)
            model.base = base/base.sum()
            for depth_slope in (False, True):
                gamma = .2*np.sin(np.arange(40)/4) if depth_slope else np.zeros(40)
                model.gamma = gamma-model.base@gamma
                for offset in (np.zeros(40), np.full(40, 7.), np.r_[np.zeros(20), np.full(20, np.log(2))]):
                    model.offset = offset.copy()
                    for mean_depth, variance in ((0., 0.), (.25, 0.), (-.25, .3)):
                        moments = Lk[:, None]*np.array([[mean_depth, variance+mean_depth**2]])
                        before_base, before_gamma = model.base.copy(), model.gamma.copy()
                        expected = independent_expected(model.base, offset, model.gamma, Lk, moments)
                        old = model.profiles(counts, Lk, moments)
                        with hook.conditioned_profiles([], {}):
                            actual = model.profiles(counts, Lk, moments)
                        error = float(np.max(np.abs(actual['E']-expected)))
                        assert error <= budget
                        assert abs(float(actual['E'].sum())-float(Lk.sum())) <= budget
                        assert np.array_equal(model.base, before_base) and np.array_equal(model.gamma, before_gamma)
                        if not np.any(offset):
                            for key in old:
                                assert np.array_equal(old[key], actual[key]), key
                        model.offset = offset+7
                        with hook.conditioned_profiles([], {}):
                            gauge = model.profiles(counts, Lk, moments)
                        assert np.max(np.abs(actual['E']-gauge['E'])) <= budget
                        cases.append(dict(baseline=baseline, depth_slope=depth_slope,
                            nonconstant_offset=bool(np.any(offset != offset[0])), depth_mean=mean_depth,
                            depth_variance=variance, max_expected_error=error, budget=budget))
        # Trigger an error inside the native profile body, after conditioning.
        model.offset = np.arange(40)/100
        before_base, before_gamma = model.base.copy(), model.gamma.copy()
        try:
            with hook.conditioned_profiles([], {}):
                model.profiles(np.zeros((1, 39)), Lk, XLk)
        except (ValueError, RuntimeError) as error:
            exception = type(error).__name__
        else:
            raise AssertionError('Exception restoration probe did not raise')
        assert np.array_equal(model.base, before_base) and np.array_equal(model.gamma, before_gamma)
        (a.out/'report.json').write_text(json.dumps(dict(all_checks_passed=True, cases=cases,
            known_offset_no_dose_legacy_expected=legacy['E'].tolist(),
            known_offset_no_dose_corrected_expected=corrected['E'].tolist(), observed_integer_mean=counts.tolist(),
            exception_restoration=exception, Decimal_precision=80, biological_specimens=0,
            production_changed=False, scientific_criteria_changed=False,
            scope='Conditional mean arithmetic and native profile execution; not DNA or clone validation'), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU native profile expectations and independent Decimal80 means')


if __name__ == '__main__':
    main()
