"""Observe actual merge populations and independently audit their pooled counts."""
import argparse
from contextlib import contextmanager
from dataclasses import asdict
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from haplodose import clones
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@contextmanager
def capture(directory, records):
    original = clones.SpecimenModel.differences
    matched_code = clones.SpecimenModel.matched_differences.__code__
    merge_code = clones.SpecimenModel.merge_unsupported.__code__

    def observe(model, Ya, La, XLa, Yb, Lb, XLb, **kwargs):
        parent = sys._getframe(1)
        merge_frame = parent.f_back
        while merge_frame is not None and merge_frame.f_code is not merge_code:
            merge_frame = merge_frame.f_back
        if parent.f_code is not matched_code or merge_frame is None:
            return original(model, Ya, La, XLa, Yb, Lb, XLb, **kwargs)
        loc, merge = parent.f_locals, merge_frame.f_locals
        result = original(model, Ya, La, XLa, Yb, Lb, XLb, **kwargs)
        assert model.alleles is None and model.graph is None and model.rows.n <= 600
        path = directory/f'comparison_{len(records):04d}.npz'
        Qa, Qb = kwargs['count_variance']
        ca, cb = kwargs['weight_scale']
        np.savez_compressed(path, pool_weights=loc['W'].toarray(),
            test=model.rows.take(np.arange(model.rows.n), model.test_fold),
            scale=model.fold_scale[model.test_fold], L=model.L, x=model.x,
            rows=model.rows.rows, labels=loc['labels'], assignment=merge['assignment_R'],
            full_assignment=merge['R'], doses=np.asarray(merge['doses']),
            Ya=Ya, La=La, XLa=XLa, Yb=Yb, Lb=Lb, XLb=XLb,
            Qa=Qa, Qb=Qb, ca=ca, cb=cb)
        source = Path(model.rows.Y.filename).parent
        seeds = [int(p.removeprefix('seed_')) for p in source.parts if p.startswith('seed_')]
        assert len(seeds) == 1
        records.append(dict(arrays=str(path), seed=seeds[0], bins=str(source/'bins.tsv'),
            pairs=loc['pairs'], results=result, params=asdict(model.diff_pp)))
        return result

    clones.SpecimenModel.differences = observe
    try:
        yield
    finally:
        clones.SpecimenModel.differences = original


def trace(args):
    runner = load('merge_population_runner', Path(__file__).with_name('membership_stage_panel.py'))
    parity = load('merge_population_parity', Path(__file__).with_name('split_draw_stability.py')).parity
    records, argv = [], sys.argv
    (args.out/'arrays').mkdir()
    try:
        sys.argv = [str(runner.__file__), '--out', str(args.out/'panel'), '--stage', 'hard',
            '--panel-root', str(args.task.parent/'recovery_20261107/recovery_panel_v1'),
            '--assertion-source', str(args.task.parent/'recovery_20261107/baseline/tests/test_split_mode.py')]
        with capture(args.out/'arrays', records):
            runner.main()
        parity(args.out/'panel', args.task.parent/'exact_profile_merge_20261109/production_v1/panel')
        assert records, 'Observer captured no actual merge comparisons'
        return dict(all_public_outputs_identical=True, all30_criteria_identical=True,
            merge_comparisons=len(records), production_changed=False, biological_specimens=0)
    finally:
        write_json(records, args.out/'comparisons.json')
        sys.argv = argv


def audit(args):
    trace_dir = args.task/args.trace_directory
    manifest = json.loads((trace_dir/'run_manifest.json').read_text())
    assert manifest['status'] == 'complete'
    for group in ('input', 'output'):
        for record in manifest[group]['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
    assert json.loads((trace_dir/'report.json').read_text())['all_public_outputs_identical']
    reports, checks = [], 0
    for item in json.loads((trace_dir/'comparisons.json').read_text()):
        z = np.load(item['arrays'])
        weights, counts = z['pool_weights'], z['test']
        pooled = np.stack([v for pair in zip(z['Ya'], z['Yb']) for v in pair])
        variance = np.stack([v for pair in zip(z['Qa'], z['Qb']) for v in pair])
        exposure = np.array([v for pair in zip(z['La'], z['Lb']) for v in pair])
        moments = np.stack([v for pair in zip(z['XLa'], z['XLb']) for v in pair])
        for group in range(weights.shape[1]):
            w = weights[:, group]
            for j in range(counts.shape[1]):
                for expected, power in ((pooled[group, j], 1), (variance[group, j], 2)):
                    observed = math.fsum(float(a)**power*float(b) for a, b in zip(w, counts[:, j]))
                    guard = 32*np.finfo(float).eps*len(w)*max(abs(observed), 1.)
                    assert abs(observed-expected) <= guard
                    checks += 1
            observed = math.fsum(float(a)*float(b) for a, b in zip(w, z['L']))*float(z['scale'])
            assert abs(observed-exposure[group]) <= 32*np.finfo(float).eps*len(w)*max(abs(observed), 1.)
            for power in (1, 2):
                terms = [float(a)*float(b)*float(c)**power*float(z['scale'])
                    for a, b, c in zip(w, z['L'], z['x'])]
                guard = 32*np.finfo(float).eps*len(w)*max(math.fsum(abs(t) for t in terms), 1.)
                assert abs(math.fsum(terms)-moments[group, power-1]) <= guard
        seed_dir = args.task.parent/f'recovery_20261107/recovery_panel_v1/seed_{item["seed"]}'
        truth = pd.read_csv(seed_dir/'data/cells.tsv', sep='\t')
        metadata = pd.read_csv(seed_dir/'model_cells.tsv', sep='\t')
        assert truth.cell_id.tolist() == metadata.cell_id.tolist()
        selected = truth.iloc[z['rows']]
        names = selected.truth_clone.to_numpy()
        bins = pd.read_csv(item['bins'], sep='\t')
        region = bins.chromosome.eq('chr12').to_numpy() & (bins.start.to_numpy() < 10000000)
        for i, (a, b) in enumerate(item['pairs']):
            populations = {}
            contrasts = {}
            for mode in ('full_MAP', 'unweighted', 'depth_matched'):
                if mode == 'full_MAP':
                    full_labels = z['full_assignment'].argmax(axis=1)
                    wa, wb = (full_labels == a).astype(float), (full_labels == b).astype(float)
                elif mode == 'unweighted':
                    wa, wb = (z['labels'] == a).astype(float), (z['labels'] == b).astype(float)
                else:
                    wa, wb = weights[:, 2*i], weights[:, 2*i+1]
                if min(wa.sum(), wb.sum()) <= 0:
                    populations[mode] = dict(untestable=True, clone_fractions=None,
                        clone1_separation=None, weight_totals=[float(wa.sum()), float(wb.sum())])
                    contrasts[mode] = None
                    continue
                fractions = [[float(w[names == label].sum()/w.sum()) for label in ('clone1', 'clone2')]
                    for w in (wa, wb)]
                populations[mode] = dict(clone_fractions=fractions,
                    clone1_separation=abs(fractions[0][0]-fractions[1][0]),
                    weight_totals=[float(wa.sum()), float(wb.sum())])
                ya, yb = wa@counts, wb@counts
                la, lb = wa@z['L']*z['scale'], wb@z['L']*z['scale']
                # Descriptive common-support pseudocount log-rate contrast; no inference claim.
                d = np.log((ya+.5)/(la+.5*counts.shape[1]))-np.log((yb+.5)/(lb+.5*counts.shape[1]))
                d -= np.median(d)
                contrasts[mode] = float(np.mean(d[region]))
            reports.append(dict(seed=item['seed'], specimen=str(selected.specimen.iloc[0]),
                comparison=item['arrays'], pair=[a,b], populations=populations,
                chr12_first10Mb_centered_mean=contrasts, returned_segments=item['results'][i]['segments']))
    assert reports and checks > 0, 'A zero-comparison audit cannot establish validity'
    write_json(reports, args.out/'comparisons.json')
    return dict(all_checks_passed=True, scalar_count_checks=checks, pairs=len(reports),
        inference_changed=False, biological_specimens=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--mode', choices=('trace', 'audit'), required=True)
    parser.add_argument('--trace-directory', default='trace_v2')
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md')
    if (args.task/'OBSERVER_AMENDMENT.md').exists():
        inputs['amendment'] = args.task/'OBSERVER_AMENDMENT.md'
    if (args.task/'MEMBERSHIP_FOLLOWUP.md').exists():
        inputs['membership_followup'] = args.task/'MEMBERSHIP_FOLLOWUP.md'
    if args.mode == 'audit':
        inputs['trace_manifest'] = args.task/args.trace_directory/'run_manifest.json'
    rec = start_run(args.out, inputs, dict(mode=args.mode, trace_directory=args.trace_directory),
        dict(scope='Conditional merge population diagnostic; no method change'), mode='merge_population')
    try:
        write_json(trace(args) if args.mode == 'trace' else audit(args), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU observed production/scalar independent pools')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
