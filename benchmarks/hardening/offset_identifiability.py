"""Shared offset and common tumour dosage require an external anchor."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.special import softmax
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--archived-helper', type=Path)
    p.add_argument('--expect-failure', action='store_true')
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    helper_path = a.archived_helper or root/'offset_nb_objective.py'
    helper = load('objective', helper_path)
    reference = load('reference', root/'offset_nb_numerics.py')
    fixture = load('fixture', root/'offset_expectation.py')
    rec = start_run(a.out, dict(script=Path(__file__), helper=helper_path,
        reference=root/'offset_nb_numerics.py', fixture=root/'offset_expectation.py',
        design=a.task/'IDENTIFICATION_DESIGN.md'), {},
        dict(scope='Conditional likelihood identification counterexample; biological n=0'), mode='offset_identifiability')
    try:
        cases = []
        for condition in ('none', 'depth', 'programme', 'both'):
            model, _, _, _, _, _ = fixture.fixture(condition, 'hard', physical_depth=True)
            model.phi, model.offset = np.array([.2, .35]), np.zeros(2)
            doses = np.array([[0., 0.], [np.log(2), 0.], [-np.log(2), 0.]])
            for normal_mass in (0., .01, .5, 1.):
                R = np.array([[normal_mass, 1-normal_mass, 0.]]*2+[[0., 0., 1.]]*2)
                data = helper.arrays(model, R, doses)
                for alpha in (-1., -.5, -.1, .1, .5, 1.):
                    offset = np.array([alpha, -alpha])
                    proposed = doses.copy()
                    proposed[1:] -= offset
                    # Exact algebra: for each tumour k>0, offset+dose'_k=dose_k.
                    eta0 = np.log(data['base'])+doses[1:, None]+data['nuisance'][None]
                    eta1 = np.log(data['base'])+offset+proposed[1:, None]+data['nuisance'][None]
                    discrepancy = float(np.max(np.abs(softmax(eta0, axis=-1)-softmax(eta1, axis=-1))))
                    forward_budget = 64*data['Y'].size*np.finfo(float).eps
                    assert discrepancy <= forward_budget
                    gain, _, budget = helper.gain_gradient(data, offset, proposed)
                    decimal = reference.decimal_gain(data, offset, proposed)
                    valid = (abs(gain-decimal) <= budget and
                             (abs(gain) <= budget if normal_mass == 0 else gain < -budget))
                    cases.append(dict(condition=condition, reference_responsibility=normal_mass,
                        alpha=alpha, gain=gain, Decimal_gain=decimal, budget=budget,
                        tumour_probability_difference=discrepancy, forward_budget=forward_budget,
                        checks_passed=bool(valid)))
        # Only reference-weighted observations carry information along this direction.
        for condition in ('none', 'depth', 'programme', 'both'):
            for alpha in (-1., -.5, -.1, .1, .5, 1.):
                rows = [v for v in cases if v['condition'] == condition and v['alpha'] == alpha]
                full = next(v for v in rows if v['reference_responsibility'] == 1.)
                for row in rows:
                    valid = abs(row['gain']-row['reference_responsibility']*full['gain']) <= row['budget']+full['budget']
                    row['mass_scaling_passed'] = bool(valid)
        passed = all(row['checks_passed'] and row['mass_scaling_passed'] for row in cases)
        (a.out/'report.json').write_text(json.dumps(dict(cases=cases, all_checks_passed=passed,
            conditional_Q_cases=96, tumour_probability_vectors=768, tumour_mean_components=1536,
            biological_specimens=0, production_changed=False, scientific_criteria_changed=False,
            conclusion='RNA common tumour dosage and free shared offset are indistinguishable when reference-component weight is zero; fitted reference responsibilities alone do not establish an independent normal anchor'), indent=2)+'\n')
        assert (not passed) if a.expect_failure else passed
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU normalized forward equivalence, Decimal80 Q and anchor-mass scaling')


if __name__ == '__main__':
    main()
