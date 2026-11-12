"""Matched-admixture null calibration of the frozen whole-workflow statistic.

Research only: one malignant CNA population plus normal cells, with or without
regional continuous RNA state. No biological calibration claim.
"""
import argparse
from copy import deepcopy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import sys
from fractions import Fraction

import jax
import numpy as np
import pandas as pd
from scipy import sparse

from haplodose import clones, pipeline
from haplodose.io import write_json
from haplodose.provenance import file_record, start_run, finish_run
from haplodose.simulate import simulate, DEFAULT_SPECIMENS


N = 86
CONDITIONS = ('matched', 'state')
SEED_BASE = 2026111031000
RANGES = ((0, 29), (29, 58), (58, 86))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def csr_hash(X):
    h = hashlib.sha256()
    for array in (X.data, X.indices, X.indptr, np.asarray(X.shape)):
        h.update(str(array.dtype).encode())
        h.update(array.tobytes())
    return h.hexdigest()


def generate(rep, condition):
    assert condition in CONDITIONS
    seed = SEED_BASE+rep
    specimens = deepcopy(DEFAULT_SPECIMENS)
    specimens[0]['clones'] = [dict(frac=1., events=[('chr3', 0., 1., 3), ('chr8', .5, 1., 1)])]
    captured = []

    def observe(frame, event, result):
        if event != 'return' or frame.f_code.co_name != 'draw':
            return
        if frame.f_code.co_filename != simulate.__code__.co_filename:
            return
        loc = frame.f_locals
        if loc['cols']['specimen'] == 'S1' and loc['cols']['truth_clone'] == 'clone1':
            captured.append((loc['L'].copy(), loc['share'].copy()))

    assert sys.getprofile() is None
    sys.setprofile(observe)
    try:
        X, genes, metadata, positions, _, _ = simulate(seed=seed, specimens=specimens)
    finally:
        sys.setprofile(None)
    assert len(captured) == 1
    exposure, share = captured[0]
    query = np.flatnonzero(metadata.specimen.eq('S1') & metadata.truth_clone.eq('clone1'))
    assert len(query) == len(exposure) == 510
    assert int((metadata.specimen.eq('S1') & metadata.truth_clone.eq('normal')).sum()) == 90
    first, stop = int(query[0]), int(query[-1])+1
    np.testing.assert_array_equal(query, np.arange(first, stop))
    preserved = sparse.vstack([X[:first], X[stop:]]).tocsr()
    base = exposure[:, None]*share
    region = positions.chromosome.eq('chr12').to_numpy() & (positions.start.to_numpy() <= 10000000)
    delta = .6*np.log(1.5)
    sigma = delta*np.sqrt(.6*.4)
    state = np.random.default_rng(np.random.SeedSequence([SEED_BASE, rep, 1])).normal(0., sigma, 510)
    effect = state if condition == 'state' else np.zeros(510)
    mean = base*np.exp(effect[:, None]*region[None])
    mean *= (exposure/mean.sum(axis=1))[:, None]
    gamma = np.random.default_rng(np.random.SeedSequence([SEED_BASE, rep, 3])).gamma(10., .1, mean.shape)
    counts = np.random.default_rng(np.random.SeedSequence([SEED_BASE, rep, 4])).poisson(mean*gamma).astype(np.int32)
    result = sparse.vstack([X[:first], sparse.csr_matrix(counts), X[stop:]]).tocsr()
    np.testing.assert_array_equal(result[:first].data, X[:first].data)
    np.testing.assert_array_equal(result[stop:].data, X[stop:].data)
    assert metadata.cell_id.is_unique and result.shape == X.shape
    arrays = dict(base_mean=base, mean=mean, exposure=exposure, state=state,
        effect=effect, region=region, query_rows=query, sigma=np.array(sigma))
    return seed, result, genes, metadata, positions, arrays, csr_hash(preserved)


def original_gate(root):
    old = root/'dev/specimen_workflow_gate_20261109'
    path = old/'source_v1/benchmarks/hardening/specimen_workflow_gate.py'
    manifest = json.loads((old/'bank_v1/run_manifest.json').read_text())
    assert manifest['status'] == 'complete'
    assert file_record(path)['sha256'] == manifest['input']['files']['script']['sha256']
    return load('admixture_unchanged_statistic', path), path


def preflight(args):
    assert jax.default_backend() == 'cpu'
    # Retain the existing86-draw maximum/rank rule; the older four-family
    # tolerance derivation is stronger than the two-family calibration here.
    assert 4*Fraction(19, 20)**N <= Fraction(1, 20)
    assert 4*Fraction(19, 20)**(N-1) > Fraction(1, 20)
    assert 2*Fraction(19, 20)**N <= Fraction(1, 20)
    gate, _ = original_gate(args.task.parent.parent)
    assert gate.decision(3., [np.full(N, 2.), np.ones(N)])['passed']
    assert not gate.decision(2., [np.full(N, 2.), np.ones(N)])['passed']
    assert not gate.decision(0., [np.zeros(N), np.zeros(N)])['passed']
    generated = [generate(1000, condition) for condition in CONDITIONS]
    assert generated[0][-1] == generated[1][-1]
    assert generated[0][3].equals(generated[1][3])
    assert generated[0][4].equals(generated[1][4])
    checks = []
    for condition, item in zip(CONDITIONS, generated):
        _, X, genes, metadata, positions, arrays, digest = item
        assert X.dtype.kind in 'iu' and np.all(X.data >= 0)
        base, mean, L, effect, region = (arrays[k] for k in ('base_mean', 'mean', 'exposure', 'effect', 'region'))
        arithmetic = 64*len(genes)*np.finfo(float).eps
        np.testing.assert_array_equal(arrays['query_rows'], generated[0][5]['query_rows'])
        assert np.isfinite(mean).all() and np.all(mean > 0)
        assert np.max(np.abs(mean.sum(axis=1)-L)/L) <= arithmetic
        # Independent scalar normalization and log contrast, at fixed rows and
        # the first regional/nonregional genes. No distributional fit or truth
        # score chooses these checks.
        inside, outside = int(np.flatnonzero(region)[0]), int(np.flatnonzero(~region)[0])
        scalar = []
        for cell in (0, 254, 509):
            denominator = math.fsum(float(b)*math.exp(float(effect[cell])) if hit else float(b)
                                   for b, hit in zip(base[cell], region))
            for gene in (inside, outside):
                expected = float(L[cell])*float(base[cell, gene])*math.exp(float(effect[cell])*bool(region[gene]))/denominator
                assert abs(float(mean[cell, gene])-expected) <= arithmetic*expected
            contrast = math.log(float(mean[cell, inside])/float(mean[cell, outside]))-math.log(float(base[cell, inside])/float(base[cell, outside]))
            assert abs(contrast-float(effect[cell])) <= arithmetic*max(1., abs(float(effect[cell])))
            scalar.append(dict(cell=cell, scalar_log_contrast=contrast, expected=float(effect[cell])))
        checks.append(dict(condition=condition, shape=X.shape, unmodified_CSR_sha256=digest,
            query_tumour_cells=510, query_normal_cells=90, relative_arithmetic_budget=arithmetic,
            independent_scalar_checks=scalar))
    write_json(dict(passed=True, conditions=checks, null_bank_size=N,
        confidence_argument='2*(19/20)^86 <= .05; fixed complete-workflow law only',
        evaluation_outcomes_read=False, biological_specimens=0), args.out/'report.json')


def worker(args):
    assert args.condition in CONDITIONS and (args.start, args.stop) in RANGES
    root = args.task.parent.parent
    check = args.task/'preflight_v1'
    assert json.loads((check/'run_manifest.json').read_text())['status'] == 'complete'
    assert json.loads((check/'report.json').read_text())['passed']
    gate, _ = original_gate(root)
    writer = load('admixture_frozen_input_writer', Path(__file__).with_name('gene_veto_fresh.py'))
    records = []
    for rep in range(args.start, args.stop):
        seed, X, genes, metadata, positions, arrays, digest = generate(rep, args.condition)
        case = args.out/f'seed_{seed}';case.mkdir()
        writer.write_input(case/'input', X, genes, metadata, positions)
        np.savez_compressed(case/'generation.npz', **arrays)
        pipeline.run(case/'input/counts.h5', case/'input/model_cells.tsv', case/'input/positions.tsv',
            case/'fit', params=clones.CloneParams(split_mode='thin'), chunk=97,
            counts_kind='raw', genome_build='hg38', coordinate_system='one-based-closed', gene_namespace='symbol')
        assert json.loads((case/'fit/run_manifest.json').read_text())['status'] == 'complete'
        value = gate.statistic(case/'fit')
        assert math.isfinite(value) and value >= 0
        records.append(dict(condition=args.condition, replicate=rep, seed=seed,
            source_unmodified_CSR_sha256=digest, case=str(case), score=value,
            fit_manifest=file_record(case/'fit/run_manifest.json'),
            input_counts=file_record(case/'input/counts.h5')))
        write_json(records, args.out/'cases.json')
        print(json.dumps(dict(condition=args.condition, fits_complete=len(records),
            total=args.stop-args.start, threshold_selected=False)), flush=True)
    write_json(dict(passed=True, cases=records, count=args.stop-args.start,
        truth_removed_from_inference=True, calibration_only=True), args.out/'report.json')


def bank(args):
    records = []
    for condition in CONDITIONS:
        for first, stop in RANGES:
            folder = args.task/f'cal_{condition}_{first:02d}_v1'
            assert json.loads((folder/'run_manifest.json').read_text())['status'] == 'complete'
            report = json.loads((folder/'report.json').read_text())
            assert report['passed'] and report['count'] == stop-first
            assert [r['replicate'] for r in report['cases']] == list(range(first, stop))
            for row in report['cases']:
                for key in ('fit_manifest', 'input_counts'):
                    assert file_record(row[key]['path'])['sha256'] == row[key]['sha256']
            records.extend(report['cases'])
    assert len(records) == 2*N
    banks = {}
    for condition in CONDITIONS:
        rows = sorted((r for r in records if r['condition'] == condition), key=lambda r:r['replicate'])
        assert [r['replicate'] for r in rows] == list(range(N))
        assert [r['seed'] for r in rows] == list(range(SEED_BASE, SEED_BASE+N))
        banks[condition] = [r['score'] for r in rows]
    by_id = {(r['condition'], r['replicate']):r for r in records}
    for rep in range(N):
        assert by_id['matched', rep]['source_unmodified_CSR_sha256'] == by_id['state', rep]['source_unmodified_CSR_sha256']
    write_json(banks, args.out/'banks.json')
    write_json(dict(passed=True, condition_datasets=2*N, independent_base_seeds=N,
        complete_workflows=2*N, thresholds={k:max(v) for k,v in banks.items()},
        conditional_calibration_family=2,
        calibration_tolerance_failure_probability_upper=float(2*Fraction(19, 20)**N),
        empirical_error_power_validation='not run', original_recovery_probe='not run',
        biological_specimens=0, production_changed=False), args.out/'report.json')


def recovery(args):
    root = args.task.parent.parent
    bank_dir = args.task/'bank_v1'
    assert json.loads((bank_dir/'run_manifest.json').read_text())['status'] == 'complete'
    banks = json.loads((bank_dir/'banks.json').read_text())
    banks = [banks['matched'], banks['state']]
    gate, _ = original_gate(root)
    scorer = load('admixture_original_criteria', Path(__file__).with_name('recovery_baseline.py'))
    tests = load('admixture_unchanged_assertion', root/'dev/recovery_20261107/baseline/tests/test_split_mode.py')
    baseline = root/'dev/shared_profile_merge_20261110/probe_v1/baseline'
    rows = []
    assert json.loads((baseline/'run_manifest.json').read_text())['status'] == 'complete'
    for before in json.loads((baseline/'summary.json').read_text())['cases']:
        seed = before['seed'];folder = baseline/f'seed_{seed}'
        value = gate.statistic(folder/'fit')
        decision = gate.decision(value, banks)
        ranks = [Fraction(1+sum(s >= value for s in b), N+1) for b in banks]
        assert decision['passed'] == (max(ranks) <= Fraction(1, N+1))
        assert decision['passed'] == (value > max(max(b) for b in banks))
        cells = pd.read_csv(folder/'assertion_view/cells.tsv.gz', sep='\t')
        differences = pd.read_csv(folder/'fit/differences.tsv', sep='\t')
        profiles = pd.read_csv(folder/'fit/clones.tsv', sep='\t')
        assert scorer.criteria(cells, differences) == before['criteria']
        names = set(profiles.loc[profiles.specimen.eq('S1') & ~profiles.normal_like, 'clone'])
        if len(names) > 1 and not decision['passed']:
            selected = cells.specimen.eq('S1') & cells.malignant_call
            cells.loc[selected, 'clone'] = 'S1_unresolved'
            cells.loc[selected, 'status'] = 'unresolved_workflow_evidence'
            differences = differences.loc[~(differences.specimen.eq('S1')
                & differences.clone_a.isin(names) & differences.clone_b.isin(names))]
        view = args.out/f'seed_{seed}';view.mkdir()
        cells.to_csv(view/'cells.tsv.gz', sep='\t', index=False)
        differences.to_csv(view/'differences.tsv', sep='\t', index=False)
        shutil.copyfile(folder/'fit/qc.json', view/'qc.json')
        passed = True
        try:
            tests.test_thin_records_mode_and_recovers_clones({'thin':view})
        except AssertionError:
            passed = False
        after = scorer.criteria(cells, differences)
        rows.append(dict(seed=seed, score=value, gate=decision,
            raw_passed=before['original_assertion_passed'], gated_passed=passed,
            raw_criteria=before['criteria'], gated_criteria=after,
            criterion_losses=[k for k,v in before['criteria'].items() if v['passed'] and not after[k]['passed']]))
    assert len(rows) == 5
    retained = not any(row['criterion_losses'] for row in rows)
    write_json(dict(cases=rows, raw_recovered=sum(r['raw_passed'] for r in rows),
        gated_recovered=sum(r['gated_passed'] for r in rows),
        all_previously_passing_criteria_retained=retained,
        continuation_gate_passed=retained, empirical_error_power_validation='not run',
        original_three_recovery_failures_fixed=False, biological_specimens=0,
        scientific_criteria_changed=False, production_changed=False, promotion=False), args.out/'report.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--mode', choices=('preflight', 'worker', 'bank', 'recovery'), required=True)
    parser.add_argument('--condition', choices=CONDITIONS)
    parser.add_argument('--start', type=int)
    parser.add_argument('--stop', type=int)
    args = parser.parse_args()
    _, gate_path = original_gate(args.task.parent.parent)
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md', statistic=gate_path)
    if args.mode == 'worker':
        inputs['numeric_gate'] = args.task/'preflight_v1/report.json'
        inputs['input_writer'] = Path(__file__).with_name('gene_veto_fresh.py')
    if args.mode == 'bank':
        for condition in CONDITIONS:
            for first, stop in RANGES:
                folder = args.task/f'cal_{condition}_{first:02d}_v1'
                inputs[f'{condition}{first}_report'] = folder/'report.json'
                inputs[f'{condition}{first}_manifest'] = folder/'run_manifest.json'
    if args.mode == 'recovery':
        inputs['banks'] = args.task/'bank_v1/banks.json'
        inputs['bank_manifest'] = args.task/'bank_v1/run_manifest.json'
        baseline = args.task.parent.parent/'dev/shared_profile_merge_20261110/probe_v1/baseline'
        inputs['original_summary'] = baseline/'summary.json'
        inputs['original_manifest'] = baseline/'run_manifest.json'
        inputs['scorer'] = Path(__file__).with_name('recovery_baseline.py')
        inputs['assertions'] = args.task.parent.parent/'dev/recovery_20261107/baseline/tests/test_split_mode.py'
        for seed in (20261031, 20261107, 20261108, 20261109, 20261110):
            for name in ('cells.tsv.gz', 'clones.tsv', 'differences.tsv', 'qc.json'):
                inputs[f'{seed}_{name}'] = baseline/f'seed_{seed}/fit'/name
            inputs[f'{seed}_assertion_view'] = baseline/f'seed_{seed}/assertion_view/cells.tsv.gz'
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, inputs, dict(mode=args.mode, condition=args.condition,
        start=args.start, stop=args.stop, seed_base=SEED_BASE, bank_size=N),
        dict(scope='Matched-admixture complete-workflow calibration; biological n=0'),
        mode='admixture_workflow_gate_'+args.mode)
    try:
        {'preflight':preflight, 'worker':worker, 'bank':bank, 'recovery':recovery}[args.mode](args)
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU isolated installed complete native workflow')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
