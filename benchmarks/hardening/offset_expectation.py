"""Independent rational pools for the shared zero-offset conditional residual."""
import argparse
import ast
import importlib.util
from fractions import Fraction
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from haplodose.clones import Rows, SpecimenModel
from haplodose.provenance import start_run, finish_run


def fixture(condition, membership, signal=False, physical_depth=False):
    depth = condition in ('depth', 'both')
    programme = condition in ('programme', 'both')
    L = np.array([12240., 3060.]*2) if physical_depth else np.full(4, 3060.)
    x = np.log(L)-np.log(L).mean() if physical_depth else np.array([math.log(2), -math.log(2)]*2)
    gamma = np.array([.5, -.5]) if depth else np.zeros(2)
    U = np.array([[math.log(2)], [-math.log(2)]]*2) if programme else np.zeros((4, 1))
    W = np.array([[1., -1.]])
    doses = np.array([[0., 0.], [math.log(2), 0.]])
    R = np.full((4, 2), .5) if membership == 'soft' else np.eye(2)[[0, 0, 1, 1]]
    means = []
    for i in range(4):
        nuisance = Fraction(2 if i % 2 == 0 else 1, 1 if i % 2 == 0 else 2) if depth else Fraction(1)
        if programme:
            nuisance *= Fraction(4 if i % 2 == 0 else 1, 1 if i % 2 == 0 else 4)
        profiles = []
        for k in range(2):
            odds = nuisance*(2 if k else 1)
            exposure = int(L[i])
            profiles.append([exposure*odds/(1+odds), exposure/(1+odds)])
        means.append(profiles)
    counts = [[sum(Fraction(str(R[i, k]))*means[i][k][j] for k in range(2)) for j in range(2)] for i in range(4)]
    assert all(v.denominator == 1 for row in counts for v in row)
    Y = np.array([[int(v) for v in row] for row in counts])
    if signal:
        Y += np.array([17, -17])
    keep = [0] if membership == 'dropped' else [0, 1]
    expected, observed = [], []
    for j in range(2):
        expected.append(sum(Fraction(str(R[i, k]))*means[i][k][j] for i in range(4) for k in keep))
        observed.append(sum(sum(Fraction(str(R[i, k])) for k in keep)*Y[i, j] for i in range(4)))
    raw = [math.log(float((o+Fraction(1, 2))/(e+Fraction(1, 2)))) for o, e in zip(observed, expected)]
    centre = .5*sum(sorted(raw))
    target = np.array([v-centre for v in raw])
    model = SimpleNamespace(base=np.array([.5, .5]), L=L, Ytot=Y.sum(axis=0),
        rows=Rows(Y, np.arange(4), 2), x=x, gamma=gamma, U=U, W=W, J=2, offset=np.array([3., -3.]))
    return model, R, doses[keep], keep, target, dict(expected=[str(v) for v in expected], observed=[str(v) for v in observed])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--expect-fixed', action='store_true')
    p.add_argument('--physical-depth', action='store_true')
    p.add_argument('--legacy-method', action='store_true', help='Execute the archived old method body, without editing a package')
    p.add_argument('--conditional-candidate', action='store_true', help='Test the explicitly unpromoted research helper')
    a = p.parse_args()
    if a.legacy_method and a.conditional_candidate:
        p.error('Choose one method')
    a.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=a.task/'DESIGN.md')
    if a.physical_depth:
        inputs['physical_design'] = a.task/'PHYSICAL_DEPTH_DESIGN.md'
    update = SpecimenModel._update_offset
    if a.legacy_method:
        archived = a.task/'baseline/python/haplodose/clones.py'
        inputs['archived_method_source'] = archived
        tree = ast.parse(archived.read_text())
        cls = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == 'SpecimenModel')
        node = next(node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == '_update_offset')
        namespace = dict(np=np, logsumexp=logsumexp)
        exec(compile(ast.Module(body=[node], type_ignores=[]), str(archived), 'exec'), namespace)
        update = namespace['_update_offset']
    if a.conditional_candidate:
        helper = Path(__file__).with_name('offset_conditional.py')
        inputs['conditional_candidate'] = helper
        spec = importlib.util.spec_from_file_location('conditional_offset', helper)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        update = module.update_offset
    run = start_run(a.out, inputs, dict(physical_depth=a.physical_depth, archived_method=a.legacy_method),
        dict(scope='Deterministic conditional means; biological n=0'), mode='offset_expectation')
    try:
        records = []
        for condition in ('none', 'depth', 'programme', 'both'):
            for membership in ('hard', 'dropped', 'soft'):
                for signal in (False, True):
                    model, R, doses, keep, target, exact = fixture(condition, membership, signal, a.physical_depth)
                    update(model, R, doses, keep)
                    error = float(np.max(np.abs(model.offset-target)))
                    operations = 64*model.rows.n*model.J*len(doses)
                    unit = np.finfo(float).eps
                    budget = operations*unit/(1-operations*unit)
                    if a.expect_fixed:
                        assert error <= budget
                    records.append(dict(condition=condition, membership=membership, signal=signal,
                        output=model.offset.tolist(), independent_target=target.tolist(), max_error=error,
                        numerical_budget=budget, **exact))
        (a.out/'report.json').write_text(json.dumps(dict(cases=records, biological_specimens=0,
            conditional_target_asserted=a.expect_fixed, physical_depth=a.physical_depth,
            source_mode=('archived_old_method_body' if a.legacy_method else
                         'unpromoted_conditional_helper' if a.conditional_candidate else 'installed_wheel'),
            scientific_criteria_changed=False), indent=2)+'\n')
        pd.DataFrame([{k: v for k, v in row.items() if k not in ('expected', 'observed')}
            for row in records]).to_csv(a.out/'counterexamples.tsv', sep='\t', index=False)
    except BaseException as error:
        finish_run(a.out, run, error=error)
        raise
    finish_run(a.out, run, backend='CPU installed offset update and exact rational conditional pools')


if __name__ == '__main__':
    main()
