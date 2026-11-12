"""Post hoc descriptive loss decomposition; exact rational arithmetic on frozen floats."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
from haplodose.provenance import start_run, finish_run


def exact(value):
    return Fraction.from_float(float(value))


def split_errors(errors):
    n = len(errors)
    bias = sum(errors)/n
    variance = sum((v-bias)**2 for v in errors)/n
    mse = sum(v*v for v in errors)/n
    assert mse == bias*bias+variance
    return bias, bias*bias, variance, mse


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--audit-name', default='audit_v2')
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), fit=a.task/'fit_v1/fitted.npz',
        truth=a.task/'panel_v1/targets.npz', summary=a.task/'fit_v1/report.json',
        audit=a.task/a.audit_name/'report.json')
    rec = start_run(a.out, inputs, {},
        dict(scope='Post hoc descriptive diagnostic; exact empirical identities, no tuning or population intervals'),
        mode='normal_offset_error_decomposition')
    try:
        fitted, target = np.load(inputs['fit']), np.load(inputs['truth'])
        summary = json.loads(inputs['summary'].read_text())
        assert json.loads(inputs['audit'].read_text())['all_checks_passed']
        records = []
        for row in summary['summary']:
            condition, method = row['condition'], row['method']
            errors, regional = [], []
            for rep in range(16):
                key = f'{condition}__{rep}'
                values = [exact(x)-exact(y) for x, y in zip(fitted[key+'__'+method], target[key+'__DNA_linked_RNA'])]
                errors.append(values)
                regional.append(sum(values[12:22])/10-(sum(values[:12])+sum(values[22:]))/30)
            decomposed = [split_errors([values[j] for values in errors]) for j in range(40)]
            bias2, variance, mse = (sum(v[index] for v in decomposed)/40 for index in (1, 2, 3))
            assert mse == bias2+variance
            mean, region_bias2, region_variance, region_mse = split_errors(regional)
            # Exact rational values independently resolve reduction/squaring errors
            # of the original floating computation;128 operations/cell/bin budget.
            operations = 128*16*40
            unit = np.finfo(float).eps
            budget = operations*unit/(1-operations*unit)*max(1., float(mse))
            assert abs(float(mse)-row['mean_DNA_linked_RNA_MSE']) <= budget
            if 'mean_regional_squared_error' in row:
                assert abs(float(region_mse)-row['mean_regional_squared_error']) <= budget
            records.append(dict(condition=condition, method=method, simulation_draws=16,
                per_bin_squared_mean_bias=float(bias2), per_bin_between_draw_variation=float(variance),
                per_bin_MSE=float(mse), regional_mean_error=float(mean),
                regional_squared_mean_bias=float(region_bias2), regional_between_draw_variation=float(region_variance),
                regional_MSE=float(region_mse), exact_identity=True, original_float_budget=budget))
        (a.out/'report.json').write_text(json.dumps(dict(records=records, exact_decompositions=len(records)*41,
            all_checks_passed=True, biological_specimens=0, production_changed=False,
            scientific_criteria_changed=False,
            scope='Between-draw variation includes random baseline/exposure and count variation; not an isolated sampling variance or a biological interval'), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU exact binary-rational empirical bias/variation identities')


if __name__ == '__main__':
    main()
