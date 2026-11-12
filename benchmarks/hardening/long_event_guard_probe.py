"""Observe an actual long-event fallback without changing the frozen workflow."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


SHORT = 'event_chrom12_bins210-218'
LONG = 'event_chrom12_bins210-228'


def load(name):
    spec = importlib.util.spec_from_file_location('long_probe_'+name,
        Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare(model, info, cells, item):
    arrays = load('gene_event_proposal')
    X = info['prop'][cells].toarray()
    share = info['share']
    fraction, _, support = arrays.composition(X, info['gbin'],
        info['resolved_gene_size']*(1-share), model.J)
    cols = item['cols']
    assert support[cols].all()
    D, EP, rP, genes, _ = arrays.region_arrays(X, info['gbin'], fraction,
        info['resolved_gene_size']*share, cols, item['E'])
    outside = ~np.isin(info['gbin'], cols)
    assert outside.any()
    depth = np.log(np.maximum(X[:, outside].sum(axis=1), 1.))
    parameters = load('depth_conditioned_genes').prepare(D, EP, rP, item['q'], depth,
        min_cells=model.params.min_cells)
    expected = item['E']*(model.fold_scale['test']/model.fold_scale['prop'])
    _, EC, rC, test_genes, _ = arrays.region_arrays(X, info['gbin'], fraction,
        info['resolved_gene_size']*(1-share), cols, expected)
    np.testing.assert_array_equal(genes, test_genes)
    selected = np.asarray(parameters['selected'], int)
    offsets = np.exp(np.asarray(parameters['offsets'])[selected, None])
    EC = EC[selected]*offsets
    rC = np.broadcast_to(rC, D.shape)[selected]
    data = dict(D=D[selected], EC=EC, rC=rC, ED=EC*share/(1-share),
        rD=rC*share/(1-share), strata=np.asarray(parameters['strata'], int)[selected],
        q=np.asarray(item['q'])[selected], selected=selected, genes=genes,
        global_rows=model.rows.rows[cells][selected])
    return data, parameters


@contextmanager
def observe(registry, method_path, old_short, out, records):
    if sys.getprofile() is not None:
        raise RuntimeError('Long-event observation requires a free profile hook')
    seen = set()

    def profile(frame, event, returned):
        if (event != 'return' or frame.f_code.co_name != 'scored'
                or frame.f_code.co_filename != str(method_path)):
            return
        loc = frame.f_locals
        source = loc['source']
        if source not in (SHORT, LONG) or source in seen:
            return
        model, cells, item = loc['self'], np.asarray(loc['cells'], int), loc['item']
        info = registry[id(model.rows)]
        assert info['case_seed'] == 202611098002 and len(cells) == 600
        data, parameters = prepare(model, info, cells, item)
        seen.add(source)
        if source == SHORT:
            with np.load(old_short) as z:
                assert set(z.files) == set(data)
                for key in z.files:
                    np.testing.assert_array_equal(data[key], z[key])
            records.append(dict(source=source, short_preparation_exact=True,
                                arrays=len(data), inference_changed=False))
            write_json(records, out/'records.json')
            return
        assert SHORT in seen
        assert loc['entry'].get('fallback') == 'declared_exact_state_cap'
        assert 3**len(item['cols']) > 1000000
        folder = out/'long'; folder.mkdir()
        np.savez_compressed(folder/'discovery.npz', **data)
        write_json(parameters, folder/'parameters.json')
        record = dict(source=source, bins=len(item['cols']),
            paths_if_enumerated=3**len(item['cols']), identified=parameters['identified'],
            contributing_groups=parameters['contributing_groups'],
            original_fallback_accepted=returned is not None, inference_changed=False)
        if not parameters['identified']:
            record.update(guard_passed=False, reason='insufficient_discovery_common_support')
            records.append(record); write_json(records, out/'records.json')
            return
        depth, pooled = load('depth_state_predictive'), load('profiled_prevalence')
        models = []
        for order in (24, 48, 96):
            print(json.dumps(dict(source=source, long_probe_order=order,
                                  status='discovery_fit')), flush=True)
            fixed = depth.fit(data['D'], data['ED'], data['rD'], data['strata'], data['q'], order)
            models.append(pooled.fit(data['D'], data['ED'], data['rD'], data['strata'],
                                     data['q'], fixed, shared=True))
            write_json(dict(models=models, confirmation_counts_opened=False,
                new_model_or_rule=False), folder/'DISCOVERY_MODELS_LOCKED.json')
        # Only this point opens confirmation observations for the added model.
        # Existing full-count nuisance/reference fits remain conditional inputs.
        full = info['X'][cells][:, data['genes']].toarray()[data['selected']]
        C = full-data['D']
        assert np.all(C >= 0) and np.all(C == np.floor(C))
        np.testing.assert_array_equal(data['D']+C, full)
        scores = []
        for fitted in models:
            value = pooled.predict(data['D'], C, data['ED'], data['EC'], data['rD'],
                data['rC'], data['strata'], fitted)
            arrays = {key: value.pop(key) for key in
                      ('per_cell', 'discovery_membership', 'joint_membership')}
            np.savez_compressed(folder/f"prediction_{fitted['order']}.npz", **arrays)
            scores.append(dict(order=fitted['order'], **value))
        diagnostic = 2*max(abs(scores[i]['score']-scores[i-1]['score']) for i in (1, 2))+sum(
            r['rounding_budget'] for r in scores)
        resolved = bool(abs(scores[-1]['score']) > diagnostic)
        record.update(scores=scores, numerical_uncertainty_diagnostic=diagnostic,
            resolved=resolved, guard_passed=bool(resolved and scores[-1]['score'] > 0),
            confirmation_count_conservation=True, model_fit_reads_confirmation=False)
        records.append(record)
        write_json(records, out/'records.json')

    sys.setprofile(profile)
    try:
        yield
    finally:
        sys.setprofile(None)


def parity(before, after):
    counts = dict(arrays=0, tables=0, folds=0)
    eligible = lambda p: p.suffix == '.npy' or p.name.endswith(('.tsv', '.tsv.gz'))
    for arm in ('production', 'pooled_prevalence'):
        a, b = before/arm, after/arm
        names = {p.name for p in (a/'fit').iterdir() if eligible(p)}
        assert names == {p.name for p in (b/'fit').iterdir() if eligible(p)}
        for name in sorted(names):
            if name.endswith('.npy'):
                np.testing.assert_array_equal(np.load(a/'fit'/name), np.load(b/'fit'/name))
                counts['arrays'] += 1
            else:
                pd.testing.assert_frame_equal(pd.read_csv(a/'fit'/name, sep='\t'),
                    pd.read_csv(b/'fit'/name, sep='\t'), check_exact=True)
                counts['tables'] += 1
        folds = {p.name for p in (a/'folds').glob('*.npz')}
        assert folds == {p.name for p in (b/'folds').glob('*.npz')}
        for name in sorted(folds):
            with np.load(a/'folds'/name) as x, np.load(b/'folds'/name) as y:
                assert x.files == y.files
                for key in x.files:
                    np.testing.assert_array_equal(x[key], y[key])
            counts['folds'] += 1
    assert counts == dict(arrays=4, tables=28, folds=1)
    return dict(counts, exact=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    prior = args.task.parent/'phase8_rep2_trace_20261110'
    assert json.loads((prior/'summary_v1/report.json').read_text())['exact_replay_gate_passed']
    args.out.mkdir(exist_ok=False)
    old = prior/'trace_v1/replay/seed_202611098002'
    short = old/'pooled_prevalence/guards/proposal_00000/discovery.npz'
    run = start_run(args.out, dict(script=Path(__file__), design=args.task/'DESIGN.md',
        exact_replay=prior/'summary_v1/report.json', original_results=prior/'trace_v1/replay/results.json',
        old_short=short, model=Path(__file__).with_name('profiled_prevalence.py'),
        depth=Path(__file__).with_name('depth_state_predictive.py'),
        conditioning=Path(__file__).with_name('depth_conditioned_genes.py')),
        dict(seed=202611098002, source=LONG, inference_changed=False),
        dict(scope='Exposed actual long-fallback mechanism; no independent error/power or biological claim'),
        mode='long_event_guard_probe')
    pilot = load('pooled_prevalence_pilot')
    guard = pilot.load('pooled_prevalence_workflow')
    original = guard.safeguard
    records, argv = [], sys.argv

    @contextmanager
    def context(method, registry, *a, **kw):
        with original(method, registry, *a, **kw):
            with observe(registry, Path(method.__file__), short, args.out, records):
                yield

    guard.safeguard = context
    try:
        sys.argv = [str(pilot.__file__), '--task', str(args.task.parent/'pooled_prevalence_pilot_20261110'),
            '--out', str(args.out/'replay'), '--condition', 'two_clone', '--start', '2', '--stop', '3']
        pilot.main()
        exact = parity(old, args.out/'replay/seed_202611098002')
        assert len(records) == 2 and records[0]['short_preparation_exact']
        long = records[1]
        write_json(dict(records=records, exact_output_parity=exact,
            long_guard_passed=long['guard_passed'], inference_changed=False,
            scientific_criteria_changed=False, biological_specimens=0, independent_evaluation=False,
            decision='Evidence for testing a scalable long-event confirmation path;not a recovery correction' if
                long['guard_passed'] else 'Reject missing long-event guard as a sufficient recovery mechanism'),
            args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU exact old workflow with observation-only long-event model fit')
    except BaseException as error:
        write_json(records, args.out/'records_failed.json')
        finish_run(args.out, run, error=error)
        raise
    finally:
        guard.safeguard, sys.argv = original, argv


if __name__ == '__main__':
    main()
