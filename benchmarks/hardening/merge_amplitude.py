"""Capture executed merge contrasts and test the fixed regional-effect prior."""
import argparse
from contextlib import contextmanager
from dataclasses import asdict
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from haplodose import clones
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def load(name):
    import importlib.util
    path = Path(__file__).with_name(name+'.py')
    spec = importlib.util.spec_from_file_location('merge_amplitude_'+name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@contextmanager
def capture(directory, records):
    original = clones.SpecimenModel.differences
    merge_code = clones.SpecimenModel.merge_unsupported.__code__
    difference_code = clones.difference_segments.__code__

    def observe(model, Ya, La, XLa, Yb, Lb, XLb, **kwargs):
        ancestor = sys._getframe(1)
        while ancestor is not None and ancestor.f_code is not merge_code:
            ancestor = ancestor.f_back
        if ancestor is None:
            return original(model, Ya, La, XLa, Yb, Lb, XLb, **kwargs)
        saved = []

        def profile(frame, event, result):
            if frame.f_code is difference_code and event == 'return' and result is not None:
                saved.append({key: np.asarray(frame.f_locals[key]).copy() for key in
                    ('idx', 'sub_reset', 'shift', 'd', 'v', 'tau', 'var', 'emit', 'trans', 'init', 'post', 'path')})

        assert sys.getprofile() is None
        sys.setprofile(profile)
        try:
            result = original(model, Ya, La, XLa, Yb, Lb, XLb, **kwargs)
        finally:
            sys.setprofile(None)
        assert len(saved) == 1 and len(Ya) == 1
        arrays = saved[0]
        arrays.update(Ya=Ya, La=La, XLa=XLa, Yb=Yb, Lb=Lb, XLb=XLb)
        path = directory/f'comparison_{len(records):04d}.npz'
        np.savez_compressed(path, **arrays)
        source = Path(model.rows.Y.filename).parent
        seeds = [int(p.removeprefix('seed_')) for p in source.parts if p.startswith('seed_')]
        assert len(seeds) == 1
        records.append(dict(arrays=str(path), seed=seeds[0], bins=str(source/'bins.tsv'),
            beta=float(model.beta), params=asdict(model.diff_pp), results=result))
        return result

    clones.SpecimenModel.differences = observe
    try:
        yield
    finally:
        clones.SpecimenModel.differences = original


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--mode', choices=('capture', 'panel'), required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    prior = args.task.parent/'merge_population_20261109'
    helper = load('merge_population')
    model = load('regional_amplitude_model')
    probe = load('regional_amplitude_probe')
    reference = load('confirmation_observation_audit')
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md', helper=Path(helper.__file__),
        model=Path(model.__file__), probe=Path(probe.__file__), reference=Path(reference.__file__),
        prior_comparisons=prior/'trace_v2/comparisons.json', prior_audit=prior/'audit_v3/report.json')
    if args.mode == 'panel':
        inputs['capture_manifest'] = args.task/'capture_v1/run_manifest.json'
    rec = start_run(args.out, inputs, dict(mode=args.mode, scale='one existing RNA difference step'),
        dict(scope='Conditional actual merge populations; no whole-workflow error/power or biology'), mode='merge_amplitude')
    try:
        assert json.loads(inputs['prior_audit'].read_text())['all_checks_passed']
        if args.mode == 'capture':
            (args.out/'payloads').mkdir()
            records = []
            try:
                with capture(args.out/'payloads', records):
                    report = helper.trace(args)
            finally:
                write_json(records, args.out/'comparisons.json')
            previous = json.loads(inputs['prior_comparisons'].read_text())
            assert len(records) == len(previous) == 4
            for current, old in zip(records, previous):
                assert current['seed'] == old['seed'] and current['results'] == old['results']
                a, b = np.load(current['arrays']), np.load(old['arrays'])
                for key in ('Ya', 'La', 'XLa', 'Yb', 'Lb', 'XLb'):
                    assert np.array_equal(a[key], b[key])
            report['actual_merge_payloads'] = 4
            report['exact_previous_pool_parity'] = True
        else:
            directory = args.task/'capture_v1'
            manifest = json.loads(inputs['capture_manifest'].read_text())
            assert manifest['status'] == 'complete'
            for group in ('input', 'output'):
                for record in manifest[group]['files'].values():
                    assert file_record(record['path'])['sha256'] == record['sha256']
            records = json.loads((directory/'comparisons.json').read_text())
            assert len(records) == 4
            rows = []
            for index, item in enumerate(records):
                z = np.load(item['arrays'])
                d, variance = z['d'][0], z['var'][0]
                step = float(np.float32(item['beta']*np.log(1.5)))
                point = model.infer(d, variance, z['sub_reset'], z['trans'], z['init'], step, 0.)
                emission = -.5*(d[:, None]-np.array([-step, 0., step]))**2/variance[:, None]
                expected, path = reference.hmm(emission, z['trans'], z['init'], z['sub_reset'])
                n = len(d)
                magnitude = 1+float(np.max(abs(emission)))+sum(abs(b['log_evidence_ratio']) for b in point['blocks'])
                u = np.finfo(float).eps/2
                guard = 128*n*u/(1-128*n*u)*magnitude
                discrepancy = float(np.max(abs(point['post']-expected)))
                assert discrepancy <= guard and np.array_equal(point['path'], path)
                fitted = model.infer(d, variance, z['sub_reset'], z['trans'], z['init'], step, step)
                bins = pd.read_csv(item['bins'], sep='\t')
                idx = z['idx']
                target = (bins.iloc[idx].chromosome.to_numpy() == 'chr12') & (bins.iloc[idx].start.to_numpy() < 10000000)
                starts = np.flatnonzero(z['sub_reset'])
                blocks = [np.arange(a, b) for a, b in zip(starts, np.r_[starts[1:], n])
                    if np.any(target[a:b])]
                assert len(blocks) == 1
                block = blocks[0]
                independent, allowance, quadrature = probe.grid_reference(d[block], variance[block],
                    z['trans'], z['init'], step, reference)
                difference = float(np.max(abs(independent-fitted['post'][block])))
                assert difference <= allowance
                np.savez_compressed(args.out/f'comparison_{index:04d}.npz',
                    post=fitted['post'], path=fitted['path'], point_post=point['post'], point_path=point['path'])
                rows.append(dict(seed=item['seed'], comparison=index, point_limit_error=discrepancy,
                    point_limit_guard=guard, whole_chr12_quadrature_error=difference,
                    whole_chr12_quadrature_guard=float(allowance), quadrature=quadrature,
                    native_max_first_half=float(np.max(z['post'][0, target, 0]+z['post'][0, target, 2])),
                    shared_max_first_half=float(np.max(fitted['post'][target, 0]+fitted['post'][target, 2])),
                    native_cores=reference.cores(z['post'][0], z['path'][0], z['sub_reset'], idx, bins, item['params']),
                    shared_cores=reference.cores(fitted['post'], fitted['path'], z['sub_reset'], idx, bins, item['params'])))
            write_json(rows, args.out/'comparisons.json')
            report = dict(all_checks_passed=True, actual_merge_comparisons=len(rows),
                full_workflow_recovery='not run', production_changed=False, biological_specimens=0)
        write_json(report, args.out/'report.json')
        finish_run(args.out, rec, backend='CPU actual merge capture / independent HMM and amplitude quadrature')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
