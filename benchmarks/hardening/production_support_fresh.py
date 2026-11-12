"""Frozen paired production calibration: fixed phase9,200 seeds per condition."""
import argparse
from fractions import Fraction
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import beta

from haplodose import clones, pipeline
from haplodose.io import write_json
from haplodose.provenance import file_record, start_run, finish_run


N = 200
ALPHA = Fraction(1, 60)
CONDITIONS = ('matched', 'state', 'two_clone')
VERSIONS = ('historical', 'corrected')


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exact_cdf(k, n, probability=Fraction(1, 20)):
    return sum(Fraction(math.comb(n, i))*probability**i*(1-probability)**(n-i)
               for i in range(k+1))


def upper(k, n, alpha):
    return 1. if k == n else float(beta.ppf(1-float(alpha), k+1, n-k))


def interval(k, n, alpha):
    return [0. if k == 0 else float(beta.ppf(float(alpha)/2, k, n-k+1)),
            1. if k == n else float(beta.ppf(1-float(alpha)/2, k+1, n-k))]


def numeric_checks():
    assert exact_cdf(0, 59) <= Fraction(1, 20) < exact_cdf(0, 58)
    allowed = max(k for k in range(N+1) if exact_cdf(k, N) <= ALPHA)
    assert exact_cdf(allowed, N) <= ALPHA < exact_cdf(allowed+1, N)
    # Independent scalar bisection of the exact rational binomial polynomial.
    comparisons = []
    for k in (0, allowed, allowed+1):
        lo, hi = Fraction(0), Fraction(1)
        for _ in range(60):
            middle = (lo+hi)/2
            if exact_cdf(k, N, middle) > ALPHA:
                lo = middle
            else:
                hi = middle
        value = upper(k, N, ALPHA)
        # beta inversion returns a float; allow an operation-count bound in the
        # inverse's probability units,not any change of the scientific threshold.
        arithmetic = (128*N+256)*np.finfo(float).eps
        assert float(lo)-arithmetic <= value <= float(hi)+arithmetic
        assert (value <= .05) == (exact_cdf(k, N) <= ALPHA)
        comparisons.append(dict(errors=k, scipy_upper=value,
            rational_root_bracket=[float(lo), float(hi)], arithmetic_budget=arithmetic))
    return dict(passed=True, N=N, phase=9, primary_family=3, alpha_per_endpoint=str(ALPHA),
        maximum_allowed_errors=allowed, zero_error59_independently_verified=True,
        exact_rational_beta_checks=comparisons)


def environment_checks(task):
    """Verify both immutable installed versions before generating any outcomes."""
    specs = json.loads((task/'ENVIRONMENTS.json').read_text())
    program = r'''
import fcntl,hashlib,importlib.metadata,json,os,sys
from pathlib import Path
import haplodose
import haplodose._core as core
import jax
s=json.loads(sys.argv[1]);prefix=Path(s['environment'])
with (prefix/'.hardening-validate.lock').open() as lock:
    fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
    expected=json.loads((Path(s['source'])/'source_manifest.json').read_text())
    package=Path(haplodose.__file__).parent
    assert package.is_relative_to(prefix) and Path(core.__file__).is_relative_to(prefix)
    for relative,digest in expected.items():
        if relative.startswith('python/haplodose/') and relative.endswith(('.py','.tsv')):
            path=package/relative.removeprefix('python/haplodose/')
            assert hashlib.sha256(path.read_bytes()).hexdigest()==digest,relative
    artifact=Path(s['artifact']);validated=json.loads(artifact.read_text())
    assert hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest()==validated['native']['sha256']
    assert hashlib.sha256(Path(validated['wheel']['path']).read_bytes()).hexdigest()==validated['wheel']['sha256']
    text=(artifact.parent/'pytest.txt').read_text()
    assert str(s['python_tests'])+' passed' in text and 'failed' not in text.strip().splitlines()[-1]
    assert 'test result: ok' in (artifact.parent/'cargo_test.txt').read_text()
    assert jax.default_backend()=='cpu'
    print(json.dumps(dict(package=str(package),wheel=validated['wheel'],native=validated['native'],
        backend=jax.default_backend(),test_summary=text.strip().splitlines()[-1],
        versions={d.metadata.get('Name',''):d.version for d in importlib.metadata.distributions()})))
'''
    runtime = Path(os.environ['HAPLODOSE_GLIBC_RUNTIME'])
    records = {}
    for version in VERSIONS:
        spec = specs[version]
        prefix = Path(spec['environment'])
        result = subprocess.run([str(runtime/'ld-linux-x86-64.so.2'), '--library-path',
            str(runtime)+':'+str(prefix/'lib'), str(prefix/'bin/python'), '-I', '-c',
            program, json.dumps(spec)], capture_output=True, text=True, check=True)
        records[version] = dict(result=json.loads(result.stdout), stderr=result.stderr)
    return records


def paths(args):
    generated = args.task/f'generated_{args.condition}_{args.start:03d}_v1'
    return generated


def preflight(task):
    report = json.loads((task/'numerics_v1/report.json').read_text())
    assert report['passed'] and report['N'] == N and report['phase'] == 9
    assert json.loads((task/'numerics_v1/run_manifest.json').read_text())['status'] == 'complete'


def generate(args, run):
    assert args.condition in CONDITIONS and args.start in (0, 100)
    preflight(args.task)
    frozen = json.loads((args.task/'GENERATOR.json').read_text())
    generator = load('production_fresh_frozen_generator', Path(frozen['path']))
    assert file_record(Path(frozen['path']))['sha256'] == frozen['sha256']
    records = []
    for rep in range(args.start, args.start+100):
        seed, X, genes, metadata, positions, arrays = generator.generate(rep, 9, args.condition)
        assert seed == 202611099000+rep
        case = args.out/f'seed_{seed}';case.mkdir()
        generator.write_input(case/'input', X, genes, metadata, positions)
        np.savez_compressed(case/'generation.npz', **arrays)
        records.append(dict(replicate=rep, seed=seed, case=str(case),
            inputs={name:file_record(case/'input'/name) for name in (
                'counts.h5', 'model_cells.tsv', 'truth_cells.tsv', 'positions.tsv')},
            generation=file_record(case/'generation.npz')))
        print(json.dumps(dict(condition=args.condition, generated=len(records))), flush=True)
    write_json(dict(condition=args.condition, phase=9, cases=records, count=100,
        generated_once=True, both_methods_same_inputs=True), args.out/'INPUTS_LOCKED.json')
    write_json(dict(passed=True, count=100, truth_removed_from_inference=True), args.out/'report.json')


def worker(args, run):
    assert args.condition in CONDITIONS and args.version in VERSIONS and args.start in (0, 100)
    preflight(args.task)
    generated = paths(args)
    assert json.loads((generated/'run_manifest.json').read_text())['status'] == 'complete'
    lock = json.loads((generated/'INPUTS_LOCKED.json').read_text())
    assert lock['count'] == 100 and lock['phase'] == 9 and lock['condition'] == args.condition
    assert [r['replicate'] for r in lock['cases']] == list(range(args.start, args.start+100))
    code = clones.SpecimenModel.em.__code__
    current_records = []

    def observe(frame, event, result):
        if event != 'return' or frame.f_code is not code:
            return
        if result is None:
            raise RuntimeError('An observed refinement did not return')
        model = frame.f_locals['self']
        doses, pi, R, LL = result[:4]
        mass = np.asarray(R).sum(axis=0)
        record = dict(components=len(doses), mass=mass.tolist(), minimum=model.params.min_cells,
            aligned=R.shape[1] == len(doses) == len(pi) == LL.shape[1],
            supported=bool(np.isfinite(mass).all() and np.all(mass[1:] >= model.params.min_cells)))
        current_records.append(record)
        assert record['aligned']
        if args.version == 'corrected':
            assert record['supported']

    output = []
    assert sys.getprofile() is None
    for record in lock['cases']:
        case = Path(record['case'])
        for item in record['inputs'].values():
            assert file_record(item['path'])['sha256'] == item['sha256']
        folder = args.out/f"seed_{record['seed']}"
        current_records = []
        sys.setprofile(observe)
        try:
            pipeline.run(case/'input/counts.h5', case/'input/model_cells.tsv',
                case/'input/positions.tsv', folder,
                params=clones.CloneParams(split_mode='thin'), chunk=97, counts_kind='raw',
                genome_build='hg38', coordinate_system='one-based-closed', gene_namespace='symbol')
        finally:
            sys.setprofile(None)
        assert current_records
        assert json.loads((folder/'run_manifest.json').read_text())['status'] == 'complete'
        output.append(dict(replicate=record['replicate'], seed=record['seed'], fit=str(folder),
            source_case=record['case'], fit_manifest=file_record(folder/'run_manifest.json'),
            returned_support=current_records))
        print(json.dumps(dict(condition=args.condition, version=args.version,
            fits_complete=len(output), truth_scored=False)), flush=True)
    write_json(dict(condition=args.condition, version=args.version, cases=output,
        count=100, truth_scored=False), args.out/'FITS_LOCKED.json')
    write_json(dict(passed=True, count=100, truth_scored=False), args.out/'report.json')


def summary(args, run):
    preflight(args.task)
    reader = load('fresh_independent_scorer', Path(__file__).with_name('combined_saved_outcomes.py'))
    assert reader.assignments([[324, 35, 1], [118, 100, 22]]) == 424
    assert reader.assignments([[360], [240]]) == 360
    scored = []
    for condition in CONDITIONS:
        for version in VERSIONS:
            for start in (0, 100):
                folder = args.task/f'{version}_{condition}_{start:03d}_v1'
                assert json.loads((folder/'run_manifest.json').read_text())['status'] == 'complete'
                lock = json.loads((folder/'FITS_LOCKED.json').read_text())
                assert lock['count'] == 100 and not lock['truth_scored']
                assert [r['replicate'] for r in lock['cases']] == list(range(start, start+100))
                for row in lock['cases']:
                    assert file_record(row['fit_manifest']['path'])['sha256'] == row['fit_manifest']['sha256']
                    truth = pd.read_csv(Path(row['source_case'])/'input/truth_cells.tsv', sep='\t')
                    value = reader.score(Path(row['fit']), truth, condition)
                    scored.append(dict(condition=condition, version=version,
                        replicate=row['replicate'], seed=row['seed'], **value))
    scores = pd.DataFrame(scored)
    assert len(scores) == 1200 and not scores.duplicated(['condition', 'version', 'replicate']).any()
    scores.to_csv(args.out/'scores.tsv', sep='\t', index=False)
    primary, descriptive = [], []
    for condition in CONDITIONS:
        for version in VERSIONS:
            subset = scores[scores.condition.eq(condition) & scores.version.eq(version)]
            assert len(subset) == N and set(subset.replicate) == set(range(N))
            candidate = version == 'corrected'
            alpha = ALPHA if candidate else Fraction(1, 20)
            if condition != 'two_clone':
                k = int(subset.null_error.sum())
                value = dict(condition=condition, version=version, errors=k, n=N,
                    upper=upper(k, N, alpha), alpha=str(alpha), endpoint='specimen_false_split',
                    exact_error_gate=bool(exact_cdf(k, N) <= alpha))
            else:
                k = int(subset.exact_chr12_recovery.sum())
                value = dict(condition=condition, version=version, recovered=k, n=N,
                    interval=interval(k, N, alpha), alpha=str(alpha), endpoint='exact_chr12_recovery')
            value.update(assignment_correct=int(subset.assignment_correct.sum()), assignment_n=600*N,
                tumour_calls=int(subset.tumour_calls.sum()), primary_family=candidate)
            (primary if candidate else descriptive).append(value)
    paired = []
    for condition in CONDITIONS:
        selected = scores[scores.condition.eq(condition)]
        a = selected[selected.version.eq('historical')].set_index('replicate').sort_index()
        b = selected[selected.version.eq('corrected')].set_index('replicate').sort_index()
        np.testing.assert_array_equal(a.seed, b.seed)
        paired.append(dict(condition=condition,
            recovery_gained=int((~a.exact_chr12_recovery & b.exact_chr12_recovery).sum()),
            recovery_lost=int((a.exact_chr12_recovery & ~b.exact_chr12_recovery).sum()),
            null_events_added=int((~a.null_error & b.null_error).sum()),
            null_events_removed=int((a.null_error & ~b.null_error).sum()),
            assignment_correct_difference=int((b.assignment_correct-a.assignment_correct).sum())))
    error_passed = all(r['exact_error_gate'] for r in primary if r['endpoint'] == 'specimen_false_split')
    write_json(dict(primary=primary, descriptive_baseline=descriptive, paired=paired,
        candidate_conditional_error_contract_passed=error_passed, primary_family=3,
        family_coverage_lower=.95, simulations_per_condition=N, condition_datasets=600,
        independent_base_seeds=200, fitted_workflows=1200, one_final_look=True,
        independent_biological_evaluation=False, biological_specimens=0,
        biological_false_split_control_established=False, genomic_membership_validated=False,
        superiority_established=False, promotion=False, scientific_criteria_changed=False),
        args.out/'report.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--mode', choices=('numerics', 'generate', 'worker', 'summary'), required=True)
    parser.add_argument('--condition', choices=CONDITIONS)
    parser.add_argument('--version', choices=VERSIONS)
    parser.add_argument('--start', type=int)
    args = parser.parse_args()
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        generator_spec=args.task/'GENERATOR.json', environments=args.task/'ENVIRONMENTS.json')
    if args.mode == 'generate':
        inputs['generator'] = Path(json.loads(inputs['generator_spec'].read_text())['path'])
    if args.mode != 'numerics':
        inputs['numeric_gate'] = args.task/'numerics_v1/report.json'
        inputs['numeric_manifest'] = args.task/'numerics_v1/run_manifest.json'
    if args.mode == 'worker':
        inputs['input_lock'] = paths(args)/'INPUTS_LOCKED.json'
        inputs['generation_manifest'] = paths(args)/'run_manifest.json'
    if args.mode == 'summary':
        inputs['scorer'] = Path(__file__).with_name('combined_saved_outcomes.py')
        for condition in CONDITIONS:
            for version in VERSIONS:
                for start in (0, 100):
                    folder = args.task/f'{version}_{condition}_{start:03d}_v1'
                    inputs[f'{version}_{condition}_{start}_lock'] = folder/'FITS_LOCKED.json'
                    inputs[f'{version}_{condition}_{start}_manifest'] = folder/'run_manifest.json'
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, inputs, dict(mode=args.mode, condition=args.condition,
        version=args.version, start=args.start, phase=9, N=N),
        dict(scope='Fixed fresh paired conditional simulation; no biological validation'),
        mode='production_support_fresh_'+args.mode)
    try:
        if args.mode == 'numerics':
            environments = environment_checks(args.task)
            write_json(dict(numeric_checks(), environments=environments), args.out/'report.json')
        else:
            {'generate':generate, 'worker':worker, 'summary':summary}[args.mode](args, run)
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU exact fixed production simulation and independent scoring')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
