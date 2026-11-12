"""Require exact installed/default parity in all181 original ordered cases."""
import argparse
import json
import os
from pathlib import Path

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    inputs = dict(script=Path(__file__), design=args.task/'REGRESSION_DESIGN.md',
        identities=root/'dev/reference_weighted_pseudobulk_20261109/regression_v1/predictions_locked.json')
    for start in (0, 46, 92, 138):
        inputs[f'part{start}'] = args.task/f'native_parity_{start:03d}_v2/report.json'
    args.out.mkdir(exist_ok=False)
    rec = start_run(args.out, inputs, dict(cases=181, equality='exact'),
        dict(scope='All original installed/default RNA regressions; no DNA read'), mode='gene_laplace_parity_summary')
    try:
        cases = []
        for start in (0, 46, 92, 138):
            report = json.loads(inputs[f'part{start}'].read_text())
            assert report['passed'] and not report['DNA_read']
            cases.extend(report['cases'])
        expected = json.loads(inputs['identities'].read_text())['cases']
        assert [r['case'] for r in cases] == [r['case'] for r in expected]
        assert [r['prediction_index'] for r in cases] == list(range(181))
        for row in cases:
            assert all(row[k] for k in ('profile_exact', 'decision_exact',
                'previous_profile_exact', 'previous_decision_exact'))
        write_json(dict(cases=181, passed=True, profile_exact=True, decision_exact=True,
            previous_profile_exact=True, previous_decision_exact=True, DNA_read=False,
            frozen_results_changed=False), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU exact exhaustive case coverage of installed parity workers')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
