"""Independent production pool arithmetic and unchanged recovery recount."""
import argparse
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from scipy import sparse
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    source = a.task/'diagnostic_v1'
    observed = json.loads((source/'report.json').read_text())
    root = Path(__file__).parent
    arithmetic = load('production_pool_arithmetic', root/'gene_veto_weights_audit.py')
    recovery = load('production_pool_recovery', root/'offset_recovery_audit.py')
    rec = start_run(a.out, dict(script=Path(__file__), design=a.task/'DESIGN.md',
        trace=source/'report.json', arithmetic=root/'gene_veto_weights_audit.py',
        recovery=root/'offset_recovery_audit.py'), {},
        dict(scope='Unchanged production; conditional pool diagnostic, biological n=0'), mode='production_pool_audit')
    argv = sys.argv
    try:
        assert observed['all_public_tables_and_posteriors_identical'] and observed['all30_criteria_identical']
        sys.argv = [str(recovery.__file__), '--out', str(a.out/'recovery'), '--task', str(a.task),
                    '--arms', 'diagnostic_v1/replay/panel']
        recovery.main()
        models, calls = {}, []
        for call in observed['calls']:
            if call['model'] not in models:
                with np.load(call['model']) as z:
                    X = sparse.csr_matrix((z['data'], z['indices'], z['indptr']), shape=tuple(z['shape'])).tocsc()
                    models[call['model']] = (X, z['sub'], z['gbin'], z['mask'], int(z['seed']))
            X, sub, gbin, mask, seed = models[call['model']]
            with np.load(call['arrays']) as z:
                labels = z['labels'][sub]
                early = bool(min(sum(labels==call['a']), sum(labels==call['b']))<10)
                assert early==call['early_cell_abstention']
                hard = [arithmetic.sums(X, (labels==k).astype(float)) for k in (call['a'], call['b'])]
                for value, key in zip(hard, ('hard_a', 'hard_b')):
                    np.testing.assert_array_equal(value, z[key])
                check = arithmetic.exact_sign(*hard, gbin, mask, call['seg'], call['min_genes'], call['alpha'])
                assert early or check['passed']==call['original_passed']
                assert call['applied_passed']==call['original_passed']
                result = dict(seed=seed, index=call['index'], caller=call['caller'],
                              original_passed=call['original_passed'], available=call['weighted']['available'],
                              early_cell_abstention=early)
                if result['available']:
                    W = z['weights'][sub]
                    moments = []
                    for j in (0, 1):
                        w = [Fraction(float(v)) for v in W[:, j]]
                        moments.append(float(sum(v*v for v in w)/sum(w)))
                    allowance = 16*(X.nnz+X.shape[0])*np.finfo(float).eps
                    np.testing.assert_allclose(z['weight_scale'], moments, rtol=allowance, atol=0)
                    effective = []
                    for j, key in enumerate(('effective_a', 'effective_b')):
                        value = arithmetic.sums(X, W[:, j])/moments[j]
                        np.testing.assert_allclose(z[key], value, rtol=allowance, atol=allowance)
                        effective.append(value)
                    check = arithmetic.exact_sign(*effective, gbin, mask, call['seg'], call['min_genes'], call['alpha'])
                    assert check['passed']==call['weighted']['passed']
                    result.update(hypothetical_passed=check['passed'], exact=check,
                                  decision_changes=not early and check['passed']!=call['original_passed'])
                else:
                    result['reason'] = call['weighted']['reason']
                calls.append(result)
        write_json(dict(all_checks_passed=True, calls=calls,
            decision_changes=[r for r in calls if r.get('decision_changes')],
            matched_available=sum(r['available'] for r in calls), total_calls=len(calls),
            original_criteria_checked=30, production_changed=False, biological_specimens=0,
            interpretation='A hypothetical aligned gene veto is not calibrated or automatically justified as a default'), a.out/'report.json')
        finish_run(a.out, rec, backend='CPU scalar sparse sums, rational moments, integer sign tails and recovery recount')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finally:
        sys.argv = argv


if __name__ == '__main__':
    main()
