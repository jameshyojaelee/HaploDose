"""Check proposal-fold mixture prevalence without changing a frozen workflow.

This is a conditional mechanism experiment, not a calibrated clone method.
Only two tumour components are considered, retaining the fitted normal mass.
"""
import argparse
from contextlib import contextmanager
from functools import lru_cache
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq
from haplodose import clones
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


@lru_cache(None)
def load(name):
    spec = importlib.util.spec_from_file_location('proposal_prevalence_'+name,
        Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fit_prevalence(LL, log_pi):
    """Exact concave scalar profile, conditional on component rates/normal mass."""
    LL, log_pi = np.asarray(LL, float), np.asarray(log_pi, float)
    if (LL.ndim != 2 or LL.shape[1] != 3 or log_pi.shape != (3,)
            or not np.isfinite(LL).all() or not np.isfinite(log_pi).all()):
        raise ValueError('Finite identified normal plus two-component likelihoods required')
    normal = math.exp(float(log_pi[0]))
    if not 0 < normal < 1:
        raise ValueError('An interior fixed normal mass is required')
    tumour_log = math.log1p(-normal)
    left = np.logaddexp(log_pi[0]+LL[:, 0], tumour_log+LL[:, 1])
    right = np.logaddexp(log_pi[0]+LL[:, 0], tumour_log+LL[:, 2])
    proportion, record = load('profiled_prevalence').prevalence(left, right)
    if not record['identified']:
        return log_pi.copy(), dict(record, changed=False)
    weights = np.array([normal, (1-normal)*(1-proportion), (1-normal)*proportion])
    with np.errstate(divide='ignore'):
        updated = np.log(weights)
    return updated, dict(record, proportion=proportion, weights=weights.tolist(),
                        fixed_normal_mass=normal, changed=True)


def numerical_checks():
    # Exact rational likelihood: two(.65-.4p) factors and one(.25+.4p).
    # The maximizer is p=1/8 with the normal probability fixed at1/5.
    LL = np.log([[.25, .75, .25], [.25, .75, .25], [.25, .25, .75]])
    old = np.log([.2, .4, .4])
    fitted, record = fit_prevalence(LL, old)
    p = record['proportion']
    eps = np.finfo(float).eps
    operations = 128*LL.size
    allowance = operations*eps/(1-operations*eps)
    error = abs(p-1/8)
    assert error <= allowance
    # Independently solve the derivative of the original probability product.
    derivative = lambda x: -.8/(.65-.4*x)+.4/(.25+.4*x)
    root = brentq(derivative, 0., 1., xtol=np.nextafter(0., 1.), rtol=4*eps)
    assert abs(p-root) <= allowance
    _, reversed_record = fit_prevalence(LL[:, [0, 2, 1]], old[[0, 2, 1]])
    assert abs(reversed_record['proportion']-(1-p)) <= allowance
    flat = LL.copy(); flat[:, 2] = flat[:, 1]
    unchanged, unidentified = fit_prevalence(flat, old)
    np.testing.assert_array_equal(unchanged, old)
    assert not unidentified['identified']
    boundary = np.log([[.25, .75, .25]])
    assert fit_prevalence(boundary, old)[1]['proportion'] == 0.
    assert fit_prevalence(boundary[:, [0, 2, 1]], old)[1]['proportion'] == 1.
    # Confirmation counts are deliberately absent from this function's inputs.
    # Each cell's additive likelihood constant cannot change the conditional fit.
    shifted, shifted_record = fit_prevalence(LL+np.array([2., -4., 8.])[:, None], old)
    assert abs(shifted_record['proportion']-p) <= allowance
    assert np.max(abs(shifted-fitted)) <= allowance
    return dict(exact_rational_maximizer='1/8', exact_error=error,
                independent_brent_error=abs(p-root), arithmetic_allowance=allowance,
                component_swap=True, flat_abstention=True, boundaries=True,
                row_constant_invariance=True, numerical_tolerance_not_scientific_threshold=True)


@contextmanager
def observe(records, directory):
    original = clones.posteriors
    merge_code = clones.SpecimenModel.merge_unsupported.__code__

    def posterior(LL, log_pi):
        result = original(LL, log_pi)
        caller = sys._getframe(1)
        if caller.f_code is not merge_code:
            return result
        model = caller.f_locals['self']
        if not model.thin or model.alleles is not None or model.graph is not None:
            raise ValueError('Only expression-only nonspatial thinning supported')
        if result.shape[1] != 3:
            records.append(dict(skipped='requires_exactly_two_tumour_components',
                                tumour_components=result.shape[1]-1))
            return result
        new_logpi, fitted = fit_prevalence(LL, log_pi)
        alternative = original(LL, new_logpi)
        entry = dict(index=len(records), fit=fitted,
                     old_weights=np.exp(log_pi).tolist(), minimum=model.params.min_cells,
                     inference_changed=False, full_count_nuisance_fixed=True,
                     rates_independently_identified=False)
        # The existing stage wrapper applies margins only when called directly
        # by merge_unsupported. Supply those SAME margins explicitly here.
        # No supported()/gene-veto call is added, hence no extra veto selection.
        for name, R in (('original', result), ('proposal_profiled', alternative)):
            labels = R.argmax(axis=1)
            counts = np.bincount(labels, minlength=3)
            row = dict(map_counts=counts.tolist(), mass=R.sum(axis=0).tolist())
            if min(counts[1:]) >= model.params.min_cells:
                pair = model.matched_differences(labels, [(1, 2)],
                    cell_weights=np.abs(R[:, 1]-R[:, 2]),
                    min_effective_cells=model.params.min_cells)[0]
                row.update(effective_cells=pair['effective_cells'],
                           testable=bool(pair.get('testable', True)),
                           untestable_reason=pair.get('untestable_reason'),
                           pre_veto_segments=pair['segments'])
            else:
                row.update(testable=False, untestable_reason='insufficient_MAP_group_size')
            entry[name] = row
        file = directory/f"merge_{entry['index']:04d}.npz"
        np.savez_compressed(file, rows=model.rows.rows, proposal_loglik=LL,
                            original_log_pi=log_pi, profiled_log_pi=new_logpi,
                            original_R=result, profiled_R=alternative)
        entry['arrays'] = str(file)
        records.append(entry)
        write_json(records, directory/'records.json')
        return result

    clones.posteriors = posterior
    try:
        yield
    finally:
        clones.posteriors = original


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    prior = args.task.parent/'phase8_rep2_trace_20261110'
    assert json.loads((prior/'summary_v1/run_manifest.json').read_text())['status'] == 'complete'
    assert json.loads((prior/'summary_v1/report.json').read_text())['exact_replay_gate_passed']
    args.out.mkdir(exist_ok=False)
    directory = args.out/'prevalence'; directory.mkdir()
    run = start_run(args.out, dict(script=Path(__file__), design=args.task/'DESIGN.md',
        original_trace=prior/'trace_v1/report.json',
        original_summary=prior/'summary_v1/report.json',
        optimizer=Path(__file__).with_name('profiled_prevalence.py')),
        dict(seed=202611098002, changed_inference=False, fixed_normal_mass=True),
        dict(scope='Exposed conditional merge-membership diagnostic; no recovery/error claim'),
        mode='proposal_prevalence_merge')
    pilot = load('pooled_prevalence_pilot')
    guard = pilot.load('pooled_prevalence_workflow')
    original = guard.safeguard
    records = []

    @contextmanager
    def context(*a, **kw):
        with original(*a, **kw):
            with observe(records, directory):
                yield

    guard.safeguard = context
    argv = sys.argv
    try:
        checks = numerical_checks()
        write_json(checks, args.out/'numerics.json')
        sys.argv = [str(pilot.__file__), '--out', str(args.out/'replay'),
            '--task', str(args.task.parent/'pooled_prevalence_pilot_20261110'),
            '--condition', 'two_clone', '--start', '2', '--stop', '3']
        pilot.main()
        before = prior/'trace_v1/replay/seed_202611098002'
        after = args.out/'replay/seed_202611098002'
        parity = load('refinement_checked_replay').exact_outputs(before, after)
        measured = [r for r in records if 'fit' in r]
        assert len(measured) == 2
        old_trace = json.loads((before/'pooled_prevalence/trace.json').read_text())
        old_merge = [r for r in old_trace['stages'] if r.get('stage') == 'merge']
        assert len(old_merge) == len(measured)
        for old, current in zip(old_merge, measured):
            np.testing.assert_array_equal(old['effective_cells'], current['original']['effective_cells'])
            assert old['testable'] == current['original']['testable']
            assert old['segments'] == len(current['original']['pre_veto_segments'])
        write_json(dict(numerics=checks, exact_output_parity=parity, records=records,
            inference_changed=False, scientific_criteria_changed=False,
            biological_specimens=0, independent_evaluation=False,
            genomic_identification=False, promotion=False), args.out/'report.json')
        finish_run(args.out, run, backend='CPU exact replay with conditional scalar prevalence profiling')
    except BaseException as error:
        write_json(records, args.out/'records_failed.json')
        finish_run(args.out, run, error=error)
        raise
    finally:
        guard.safeguard, sys.argv = original, argv


if __name__ == '__main__':
    main()
