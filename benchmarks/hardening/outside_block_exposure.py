"""Test exposure estimated outside the genomic block of each event proposal."""
import argparse
from contextlib import contextmanager
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys

import jax
import numpy as np
from scipy.special import logsumexp
from scipy.stats import nbinom

from haplodose import clones
from haplodose.io import write_json
from haplodose.provenance import file_record, start_run, finish_run


def load(name):
    path = Path(__file__).with_name(name+'.py')
    spec = importlib.util.spec_from_file_location('outside_exposure_'+name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exposure(counts, rates, sizes):
    """Solve the independent-bin NB exposure score for each row by bisection."""
    y, q, r = np.broadcast_arrays(*[np.asarray(v, float) for v in (counts, rates, sizes)])
    if (y.ndim != 2 or not y.shape[1] or not all(np.isfinite(v).all() for v in (y, q, r))
            or np.any(y < 0) or np.any(q <= 0) or np.any(r <= 0)):
        raise ValueError('Finite nonnegative counts and positive rates and sizes required')
    if np.any(y.sum(axis=1) == 0):
        return None
    lo, hi = np.zeros(len(y)), np.max(y/q, axis=1)
    if not np.isfinite(hi).all():
        raise ArithmeticError('Exposure bracket exceeds floating point range')
    # Each summand is nonnegative at zero and nonpositive at max(y/q).
    # Stop only when the remaining bracket has adjacent floating endpoints.
    for _ in range(np.finfo(float).maxexp + np.finfo(float).nmant):
        mid = lo+(hi-lo)/2
        working = (mid > lo) & (mid < hi)
        if not working.any():
            break
        mean = mid[:, None]*q
        score = np.sum((y-mean)/(1+mean/r), axis=1)
        lo = np.where(working & (score > 0), mid, lo)
        hi = np.where(working & (score <= 0), mid, hi)
    else:
        raise ArithmeticError('Exposure bracket failed to resolve')
    return lo+(hi-lo)/2


def outside_expected(counts, rates, sizes, outside):
    fitted = exposure(counts[:, outside], rates[:, outside], sizes[outside])
    return None if fitted is None else fitted[:, None]*rates


def checks():
    y = np.array([[1., 3., 9.], [2., 6., 7.]])
    q = np.tile([.125, .375, .5], (2, 1))
    r = np.array([2., 4., 3.])
    outside = np.array([True, True, False])
    actual = outside_expected(y, q, r, outside)
    # Both included bins attain their observed counts at the same exposure.
    exact = np.array([Fraction(1, 1)/Fraction(1, 8), Fraction(2, 1)/Fraction(1, 8)])
    eps = np.finfo(float).eps
    operations = 32*y.size+64
    gamma = operations*eps/(1-operations*eps)
    budget = gamma*float(max(exact))
    np.testing.assert_allclose(actual[:, -1], np.asarray(exact, float)/2, rtol=0, atol=budget)
    changed = y.copy()
    changed[:, -1] *= 100
    np.testing.assert_array_equal(outside_expected(changed, q, r, outside), actual)
    perm = np.array([2, 0, 1])
    other = outside_expected(y[::-1, perm], q[::-1, perm], r[perm], outside[perm])
    np.testing.assert_allclose(other[::-1, np.argsort(perm)], actual, rtol=0, atol=budget)
    # Independent scalar PMFs on both sides of the exact interior maximum.
    for yy, theta in zip(y, exact):
        def value(t):
            means = float(t)*q[0, outside]
            return nbinom.logpmf(yy[outside], r[outside], r[outside]/(r[outside]+means)).sum()
        assert value(theta) > max(value(theta/2), value(theta*2))
    # At theta=8 the heterogeneous fixture has terms 1/2 and -2/4.
    residuals = (Fraction(1, 2), Fraction(-2, 4))
    assert sum(residuals) == 0
    np.testing.assert_allclose(exposure([[2, 1]], [[.125, .375]], [1, 1]),
                               [8], rtol=0, atol=budget)
    assert outside_expected(np.array([[0., 0., 9.]]), q[:1], r, outside) is None
    return dict(passed=True, exact_exposures=list(map(str, exact)), rounding_budget=budget,
        independent_nb_pmf_checked=True, bin_and_cell_permutation=True,
        heterogeneous_exact_root=8,
        tested_block_count_perturbation_exact=True, zero_outside_count_abstention=True,
        upstream_nuisance_fixed=True)


@contextmanager
def replace_exposure(records, directory, *, arm):
    original = clones._fit_event_mixture
    event_code = clones.SpecimenModel.event_mixtures.__code__
    last_call, current, cache = None, {}, {}

    def fit(O, E, r, step, *, maxiter=1000):
        nonlocal last_call
        frame = sys._getframe(1)
        if frame.f_code is not event_code:
            return original(O, E, r, step, maxiter=maxiter)
        loc = frame.f_locals
        model = loc['self']
        if model.graph is not None or model.alleles is not None or model.params.split_null_draws:
            raise ValueError('Expression-only nonspatial research with no selected-parent null')
        if frame is not last_call:
            last_call = frame
            cache.clear()
            counts, sizes = loc['Yp'], loc['r_bin']
            rates = np.exp(loc['A']-logsumexp(loc['A'], axis=1, keepdims=True))
            exact = exposure(counts, rates, sizes)
            path = directory/f'call_{len(records):05d}.npz'
            payload = dict(counts=counts, rates=rates, sizes=sizes, cols=loc['cols'],
                           blocks=loc['chrom'], native_expected=loc['Ep'],
                           global_rows=model.rows.rows[np.asarray(loc['cells'], int)])
            if exact is not None:
                payload['exact_all_exposure'] = exact
            np.savez_compressed(path, **payload)
            current.update(counts=counts, rates=rates, sizes=sizes, exact=exact, path=str(path))
        block, part, cols = int(loc['c']), loc['part'], loc['cols']
        if block not in cache:
            outside = loc['chrom'] != block
            fitted = exposure(current['counts'][:, outside], current['rates'][:, outside],
                              current['sizes'][outside]) if outside.any() else None
            path = directory/f'block_{len(records):05d}.npz'
            payload = dict(outside=outside)
            if fitted is not None:
                payload['outside_exposure'] = fitted
            np.savez_compressed(path, **payload)
            cache[block] = (fitted, str(path))
        fitted, path = cache[block]
        theta, estimator = None, 'native'
        if arm != 'sham' and current['exact'] is not None:
            theta, estimator = current['exact'], 'exact_all'
        if arm == 'candidate' and fitted is not None:
            theta, estimator = fitted, 'outside_block'
        tested = E if theta is None else theta[:, None]*current['rates'][:, part]
        record = dict(block=block, first_bin=int(cols[part[0]]), last_bin=int(cols[part[-1]]),
            source_counts=str(model.rows.Y.filename), cells=len(O), call_arrays=current['path'],
            block_arrays=path, estimator=estimator, outside_supported=fitted is not None)
        records.append(record)
        try:
            result = original(O, tested, r, step, maxiter=maxiter)
        except RuntimeError as error:
            record['unresolved'] = str(error)
            raise
        minority = int(np.sum(result['responsibility'] > .5))
        record.update(gain=float(result['gain']), hard_minority=min(minority, len(O)-minority))
        return result

    clones._fit_event_mixture = fit
    try:
        yield
    finally:
        clones._fit_event_mixture = original
        last_call = None
        current.clear()
        cache.clear()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if jax.default_backend() != 'cpu':
        raise RuntimeError('Actual CPU backend required')
    args.out.mkdir(exist_ok=False)
    runner, helper = load('membership_stage_panel'), load('shared_profile_merge_probe')
    recovery = args.task.parent/'recovery_20261107'
    panel, assertions = recovery/'recovery_panel_v1', recovery/'baseline/tests/test_split_mode.py'
    run = start_run(args.out, dict(script=Path(__file__), design=args.task/'DESIGN.md',
        runner=Path(runner.__file__), helper=Path(helper.__file__), panel=panel/'manifest.json',
        assertions=assertions), dict(arms=['baseline', 'sham', 'exact_all', 'candidate']),
        dict(scope='Original exposed recovery cases with conditional outside-block exposure'),
        mode='outside_block_exposure')
    summaries, argv = {}, sys.argv
    try:
        write_json(checks(), args.out/'numerics.json')
        for arm in ('baseline', 'sham', 'exact_all', 'candidate'):
            directory = args.out/(arm+'_exposures')
            directory.mkdir()
            trace = []
            sys.argv = [str(runner.__file__), '--out', str(args.out/arm), '--panel-root', str(panel),
                        '--assertion-source', str(assertions), '--stage', 'hard']
            if arm == 'baseline':
                runner.main()
            else:
                try:
                    with replace_exposure(trace, directory, arm=arm):
                        runner.main()
                finally:
                    write_json(trace, args.out/(arm+'_trace.json'))
                assert trace
                if arm == 'candidate':
                    assert any(row['estimator'] == 'outside_block' for row in trace)
            summaries[arm] = json.loads((args.out/arm/'summary.json').read_text())['cases']
            if arm == 'sham':
                exact = helper.exact_outputs(args.out/'baseline', args.out/'sham')
                assert len(exact) == 5
                write_json(dict(passed=True, cases=exact), args.out/'SHAM_LOCKED.json')
        comparisons = []
        for before, exact, after in zip(summaries['baseline'], summaries['exact_all'],
                                        summaries['candidate'], strict=True):
            assert before['seed'] == exact['seed'] == after['seed']
            comparisons.append(dict(seed=before['seed'], before=before, exact_all=exact, after=after,
                recovered=not before['original_assertion_passed'] and after['original_assertion_passed'],
                lost=before['original_assertion_passed'] and not after['original_assertion_passed'],
                criterion_losses=[key for key, value in before['criteria'].items()
                                  if value['passed'] and not after['criteria'][key]['passed']],
                added_s3_groups=after['s3_groups'] > before['s3_groups']))
        assert len(comparisons) == 5
        full = all(r['original_assertion_passed'] for r in summaries['candidate'])
        advance = full and not any(r['lost'] or r['criterion_losses'] or r['added_s3_groups']
                                   for r in comparisons)
        write_json(dict(comparisons=comparisons,
            baseline_recovered=sum(r['original_assertion_passed'] for r in summaries['baseline']),
            candidate_recovered=sum(r['original_assertion_passed'] for r in summaries['candidate']),
            exact_all_recovered=sum(r['original_assertion_passed'] for r in summaries['exact_all']),
            continuation_passed=advance, all_original_recovered=full,
            scientific_criteria_changed=False, biological_specimens=0, production_changed=False),
            args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU installed wheel, outside-block proposal exposure')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise
    finally:
        sys.argv = argv


if __name__ == '__main__':
    main()
