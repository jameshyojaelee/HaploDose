"""Independent Decimal80 offsets and exact rational kept-population pools."""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
from haplodose.provenance import start_run, finish_run


def log_ratio(observed, expected):
    with localcontext() as context:
        context.prec = 80
        raw = []
        for o, e in zip(observed, expected):
            ratio = (o+Fraction(1, 2))/(e+Fraction(1, 2))
            raw.append((Decimal(ratio.numerator)/Decimal(ratio.denominator)).ln())
        centre = sum(raw)/2
        return [float(v-centre) for v in raw]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--before-name', default='before_v1')
    p.add_argument('--after-name', default='after_v1')
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    rec = start_run(a.out, dict(script=Path(__file__), design=a.task/'DESIGN.md',
        before=a.task/a.before_name/'report.json', after=a.task/a.after_name/'report.json'), {},
        dict(scope='Exact conditional pools, not biological calibration'), mode='offset_expectation_audit')
    try:
        before_report = json.loads((a.task/a.before_name/'report.json').read_text())
        after_report = json.loads((a.task/a.after_name/'report.json').read_text())
        assert before_report.get('physical_depth', False) == after_report.get('physical_depth', False)
        before, after = before_report['cases'], after_report['cases']
        assert len(before) == len(after) == 24
        results = []
        for old, new in zip(before, after):
            key = ('condition', 'membership', 'signal')
            assert all(old[k] == new[k] for k in key)
            condition, membership, signal = (new[k] for k in key)
            depth, programme = condition in ('depth', 'both'), condition in ('programme', 'both')
            exposure = [12240, 3060]*2 if after_report.get('physical_depth', False) else [3060]*4
            means, counts, weight = [], [], []
            for i in range(4):
                w = [Fraction(1, 2)]*2 if membership == 'soft' else [Fraction(int(i < 2)), Fraction(int(i >= 2))]
                weight.append(w)
                factor = (2 if depth else 1)*(4 if programme else 1)
                factor = Fraction(factor) if i % 2 == 0 else Fraction(1, factor)
                cells = []
                for k in range(2):
                    ratio = factor*(2 if k else 1)
                    cells.append([exposure[i]*ratio/(1+ratio), exposure[i]/(1+ratio)])
                means.append(cells)
                observed = [sum(w[k]*cells[k][j] for k in range(2)) for j in range(2)]
                if signal:
                    observed = [observed[0]+17, observed[1]-17]
                counts.append(observed)
            keep = [0] if membership == 'dropped' else [0, 1]
            obs = [sum(sum(weight[i][k] for k in keep)*counts[i][j] for i in range(4)) for j in range(2)]
            exp = [sum(weight[i][k]*means[i][k][j] for i in range(4) for k in keep) for j in range(2)]
            expected_old = [sum(sum(weight[i][k]*exposure[i] for i in range(4))*
                (Fraction(1, 2) if k == 0 else Fraction(2 if j == 0 else 1, 3)) for k in keep) for j in range(2)]
            observed_old = [sum(row[j] for row in counts) for j in range(2)]
            corrected = log_ratio(obs, exp)
            previous = log_ratio(observed_old, expected_old)
            assert all(str(v) == new['observed'][j] for j, v in enumerate(obs))
            assert all(str(v) == new['expected'][j] for j, v in enumerate(exp))
            budget = new['numerical_budget']
            assert max(abs(x-y) for x, y in zip(new['output'], corrected)) <= budget
            assert max(abs(x-y) for x, y in zip(old['output'], previous)) <= budget
            assert max(abs(x-y) for x, y in zip(new['independent_target'], corrected)) <= budget
            if not signal:
                assert corrected == [0., 0.]
            else:
                assert corrected[0] > 0 and corrected[1] < 0
            results.append(dict(condition=condition, membership=membership, signal=signal,
                before_error=max(abs(x-y) for x, y in zip(old['output'], corrected)),
                after_error=max(abs(x-y) for x, y in zip(new['output'], corrected)), budget=budget))
        (a.out/'report.json').write_text(json.dumps(dict(cases=results, Decimal_precision=80,
            exact_count_pools=48, old_new_scalar_offsets=96, all_checks_passed=True,
            biological_specimens=0, scientific_criteria_changed=False), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU Fraction pools and Decimal80 independent old/new offsets')


if __name__ == '__main__':
    main()
