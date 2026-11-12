"""Independent counts, interpolation and family decisions for combined precision."""
import argparse
from fractions import Fraction
import importlib.util
import json
from pathlib import Path

import numpy as np
from haplodose.provenance import start_run, finish_run, file_record


def verify(record):
    assert file_record(record['path'])['sha256'] == record['sha256'], record['path']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    combined = args.task / 'complete_v1'
    manifest = json.loads((combined / 'run_manifest.json').read_text())
    assert manifest['status'] == 'complete'
    for group in ('input', 'output'):
        for record in manifest[group]['files'].values():
            verify(record)
    inputs = dict(script=Path(__file__), manifest=combined / 'run_manifest.json',
                  summary=combined / 'summary.json', draws=combined / 'loss_draws.npy',
                  reference=Path(__file__).with_name('gene_factorial_precision_audit.py'))
    run = start_run(args.out, inputs, dict(family=153, draws=1000000),
                    dict(scope='Independent numerical reconstruction; biological coverage unvalidated'),
                    mode='factor_loss_summary_audit')
    try:
        spec = importlib.util.spec_from_file_location('independent_reference', inputs['reference'])
        reference = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(reference)
        summary = json.loads(inputs['summary'].read_text())
        assert summary['draws'] == 1000000 and not summary['prefix_check']
        assert summary['scientific_family'] == 153 and summary['numerical_family'] == 459
        rows = summary['comparisons']
        assert len(rows) == 153 and len({(r['endpoint'], r['a'], r['b']) for r in rows}) == 153
        draws = np.load(inputs['draws'], mmap_mode='r')
        assert draws.shape == (1000000, 51)
        probabilities = {}
        endpoints = 0
        for j, row in enumerate(r for r in rows if r['endpoint'] == 'recall_loss'):
            values = draws[:, j]
            values = values[np.isfinite(values)]
            assert len(values) == row['valid_bootstrap_replicates']
            delta = row['mean_diff']
            testable = bool(np.isfinite(delta) and len(values) >= 900000 and np.ptp(values) > 0)
            assert testable == row['testable']
            errors = values - delta
            count = sum(int(np.count_nonzero(np.abs(errors[k:k+8192]) >= abs(delta)))
                        for k in range(0, len(errors), 8192))
            probability = Fraction(count + 1, len(values) + 1) if testable else Fraction(1)
            assert row['p'] == float(probability)
            probabilities[row['endpoint'], row['a'], row['b']] = probability
            if testable:
                ordered = np.sort(errors)
                for field, q in [('ci_low', .975), ('ci_high', .025),
                                 ('family_ci_low', 1-.05/306), ('family_ci_high', .05/306)]:
                    value, budget = reference.basic_endpoint(ordered, delta, q)
                    assert abs(value - row[field]) <= budget
                    endpoints += 1
        probabilities = [probabilities.get((r['endpoint'], r['a'], r['b']),
                                          Fraction.from_float(r['p'])) for r in rows]
        # Compute Holm using sorted exact rational probabilities and a running
        # maximum, independent of the collector's pairwise implementation.
        expected = [None] * len(rows)
        running = Fraction(0)
        for rank, i in enumerate(sorted(range(len(rows)), key=probabilities.__getitem__)):
            running = max(running, min(Fraction(1), (153-rank)*probabilities[i]))
            expected[i] = running
        for row, exact in zip(rows, expected):
            budget = 3*np.finfo(float).eps*153*max(row['p'], float(exact))
            assert abs(row['p_holm'] - float(exact)) <= budget
            assert (row['p_holm'] < .05) == (exact < Fraction(1, 20))
            assert row['superior'] == bool(row['testable'] and
                                         row['p_holm'] < .05 and row['family_ci_low'] > 0)
        assert summary['resampling_rows_verified'] == 1000000
        report = dict(independent_loss_probability_counts=51,
                      independent_interval_endpoints=endpoints, exact_holm_decisions=153,
                      biological_validation=False, scientific_criteria_changed=False)
        (args.out / 'report.json').write_text(json.dumps(report, indent=2))
        for record in run['input']['files'].values():
            verify(record)
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
