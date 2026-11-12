"""Join indexed exact-loss draws and audit the unchanged 153-test family."""
import argparse
import importlib.util
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, binom
from haplodose.provenance import start_run, finish_run, file_record


def load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verified_manifest(directory):
    manifest = json.loads((directory / 'run_manifest.json').read_text())
    assert manifest['status'] == 'complete', directory
    for group in ('input', 'output'):
        for record in manifest[group]['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256'], record['path']
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--chunks', type=Path, required=True)
    parser.add_argument('--prefix-check', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    root = args.task.parents[1]
    factor = root / 'dev/gene_factorial_20261108'
    exact = factor / 'exact_uncertainty_v1'
    corr = factor / 'precision_v1'
    gate = args.task / 'prefix_audit_v1'
    total = 20000 if args.prefix_check else 1000000
    assert json.loads((gate / 'report.json').read_text())['continuation_gate_passed']
    for directory in (gate, exact, corr, factor / 'precision_audit_v1'):
        verified_manifest(directory)
    chunks = json.loads(args.chunks.read_text())
    inputs = dict(script=Path(__file__), chunks=args.chunks, design=args.task / 'DESIGN.md',
                  gate=gate / 'report.json', exact=exact / 'summary.json',
                  correlation=corr / 'summary.json',
                  correlation_audit=factor / 'precision_audit_v1/report.json',
                  scorer=Path(__file__).with_name('profile.py'),
                  reference=Path(__file__).with_name('gene_factorial_precision_audit.py'))
    cursor = 0
    manifests = []
    for i, chunk in enumerate(chunks):
        assert chunk['start'] == cursor and chunk['stop'] > cursor
        cursor = chunk['stop']
        directory = Path(chunk['path'])
        manifest = verified_manifest(directory)
        report = json.loads((directory / 'report.json').read_text())
        assert (report['start'], report['stop']) == (chunk['start'], chunk['stop'])
        manifests.append(manifest)
        for name in ('run_manifest.json', 'draws.npy', 'indices.npy', 'definition.json'):
            inputs[f'chunk{i}:{name}'] = directory / name
    assert cursor == total
    origin = manifests[0]['input']['files']
    for manifest in manifests:
        for name in ('scorer', 'cache', 'support', 'split', *[f'fit{i}' for i in range(9)]):
            assert manifest['input']['files'][name]['sha256'] == origin[name]['sha256']
    assert file_record(inputs['scorer'])['sha256'] == origin['scorer']['sha256']
    for name in ('support', 'split', *[f'fit{i}' for i in range(9)]):
        inputs[name] = Path(origin[name]['path'])
    run = start_run(args.out, inputs, dict(draws=total, scientific_family=153,
                    numerical_family=459, prefix_check=args.prefix_check),
                    dict(scope='Numerical precision only; exposed development specimens'),
                    mode='factor_loss_summary')
    try:
        base = load('profile')
        reference = load('gene_factorial_precision_audit')
        original = json.loads(inputs['exact'].read_text())['comparisons']
        rows = [dict(r) for r in original if r['endpoint'] == 'recall_loss']
        definition = json.loads(inputs['chunk0:definition.json'].read_text())
        assert len(rows) == 51 and definition['contrasts'] == [[r['a'], r['b']] for r in rows]
        draws = np.lib.format.open_memmap(args.out / 'loss_draws.npy', mode='w+',
                                          dtype=float, shape=(total, 51))
        with np.load(inputs['support'], allow_pickle=True) as data:
            loss, neutral, common, samples = (data[k] for k in ('loss', 'neutral', 'common', 'samples'))
        meta = pd.read_csv(inputs['split'], sep='\t').set_index('sample').loc[samples]
        ci = np.flatnonzero(meta.split.eq('calibration'))
        ei = np.flatnonzero(meta.split.eq('evaluation'))
        assert len(ci) == len(ei) == 24 and meta.patient.nunique() == 48
        rng = np.random.default_rng(20261106)
        scalar_checks = 0
        for j, chunk in enumerate(chunks):
            assert json.loads(inputs[f'chunk{j}:definition.json'].read_text()) == definition
            values = np.load(inputs[f'chunk{j}:draws.npy'], mmap_mode='r')
            indices = np.load(inputs[f'chunk{j}:indices.npy'], mmap_mode='r')
            count = chunk['stop'] - chunk['start']
            assert values.shape == (count, 51) and indices.shape == (count, 49)
            draws[chunk['start']:chunk['stop']] = values
            for k in range(count):
                ri = int(rng.integers(8)) + 1
                cal = rng.choice(ci, len(ci), replace=True)
                chosen = rng.choice(ei, len(ei), replace=True)
                np.testing.assert_array_equal(indices[k], np.r_[ri, cal, chosen])
                if k not in (0, count - 1):
                    continue
                with np.load(inputs[f'fit{ri}'], allow_pickle=True) as fit:
                    scored = {}
                    for method in definition['methods']:
                        threshold = base.loss_threshold(fit[method], neutral & common, cal)
                        scored[method] = np.array([r['recall_loss'] for r in
                            base.score_loss(fit[method], loss, neutral, common, threshold)])
                for column, (a, b) in enumerate(definition['contrasts']):
                    vector = (scored[a] - scored[b])[chosen]
                    vector = vector[np.isfinite(vector)]
                    if not len(vector):
                        assert np.isnan(values[k, column])
                        continue
                    expected = math.fsum(vector) / len(vector)
                    eps = np.finfo(float).eps / 2
                    gamma = (len(vector) + 1) * eps / (1 - (len(vector) + 1) * eps)
                    budget = gamma * math.fsum(abs(float(v)) for v in vector) / len(vector)
                    assert abs(expected - values[k, column]) <= budget
                    scalar_checks += 1
            print(json.dumps(dict(chunk=j, resampling_rows_verified=chunk['stop'])), flush=True)
        draws.flush()
        bounds = {}
        tail = .05 / (2 * 459)
        for j, row in enumerate(rows):
            delta = definition['point_effects'][j]
            assert delta == row['mean_diff']
            boot = np.asarray(draws[:, j])
            boot = boot[np.isfinite(boot)]
            testable = bool(np.isfinite(delta) and len(boot) >= int(.9 * total) and np.ptp(boot) > 0)
            updated = dict(row, valid_bootstrap_replicates=len(boot), testable=testable, p=1.,
                           ci_low=np.nan, ci_high=np.nan, family_ci_low=np.nan, family_ci_high=np.nan)
            bound = dict(p_lower=0., p_upper=1., numerical_status='untestable')
            if testable:
                errors = boot - delta
                count = sum(int(np.sum(np.abs(errors[k:k+10000]) >= abs(delta)))
                            for k in range(0, len(errors), 10000))
                updated['p'] = float(Fraction(count + 1, len(errors) + 1))
                ordered = np.sort(errors)
                for field, q in [('ci_low', .975), ('ci_high', .025),
                                 ('family_ci_low', 1-.05/306), ('family_ci_high', .05/306)]:
                    value = float(delta - np.quantile(errors, q))
                    independent, budget = reference.basic_endpoint(ordered, delta, q)
                    assert abs(value - independent) <= budget
                    updated[field] = value
                regions = []
                for q in (.05/306, 1-.05/306):
                    low = int(binom.ppf(tail, len(errors), q))
                    high = int(binom.isf(tail, len(errors), q)) + 1
                    assert binom.cdf(low-1, len(errors), q) <= tail
                    assert binom.sf(high-1, len(errors), q) <= tail
                    regions.append((-np.inf if low == 0 else float(ordered[low-1]),
                                    np.inf if high > len(errors) else float(ordered[high-1])))
                bound = dict(exceedances=count, n=len(errors),
                    p_lower=0. if count == 0 else float(beta.ppf(tail, count, len(errors)-count+1)),
                    p_upper=1. if count == len(errors) else float(beta.isf(tail, count+1, len(errors)-count)),
                    family_ci_low_mc=[delta-regions[1][1], delta-regions[1][0]],
                    family_ci_high_mc=[delta-regions[0][1], delta-regions[0][0]])
            if args.prefix_check:
                for field in ('p', 'testable', 'valid_bootstrap_replicates', 'ci_low', 'ci_high',
                              'family_ci_low', 'family_ci_high'):
                    assert updated[field] == row[field] or (
                        isinstance(updated[field], float) and np.isnan(updated[field]) and np.isnan(row[field]))
            rows[j] = updated
            bounds[(row['endpoint'], row['a'], row['b'])] = bound
        correlation = json.loads(inputs['correlation'].read_text())['prefixes'][-1]['comparisons']
        mapping = {(r['endpoint'], r['a'], r['b']): r for r in rows}
        whole = [mapping.get((r['endpoint'], r['a'], r['b']), dict(r)) for r in correlation]
        assert len(whole) == 153
        base.adjust_family(whole)
        independent_holm = reference.holm([Fraction.from_float(r['p']) for r in whole])
        for row, expected in zip(whole, independent_holm):
            assert abs(row['p_holm'] - float(expected)) <= 3*np.finfo(float).eps*153*max(row['p'], float(expected))
            assert (row['p_holm'] < .05) == (expected < Fraction(1, 20))
        for record in json.loads(inputs['correlation_audit'].read_text())['final_mc']:
            if record['endpoint'] != 'recall_loss':
                bounds[(record['endpoint'], record['a'], record['b'])] = record
        ordered_bounds = [dict(bounds[r['endpoint'], r['a'], r['b']],
                               endpoint=r['endpoint'], a=r['a'], b=r['b']) for r in whole]
        lower = reference.holm([r['p_lower'] for r in ordered_bounds])
        upper = reference.holm([r['p_upper'] for r in ordered_bounds])
        for record, lo, hi, row in zip(ordered_bounds, lower, upper, whole):
            record.update(holm_mc_lower=lo, holm_mc_upper=hi,
                positive_claim_numerically_stable=bool(row['superior'] and hi < .05 and
                    record.get('family_ci_low_mc', [-np.inf])[0] > 0))
        (args.out / 'summary.json').write_text(json.dumps(dict(comparisons=whole, mc=ordered_bounds,
            draws=total, scientific_family=153, numerical_family=459, prefix_check=args.prefix_check,
            independent_scalar_draw_effects=scalar_checks, resampling_rows_verified=total,
            biological_units_added=0, scope='Conditional numerical precision; biological coverage unvalidated'), indent=2))
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
