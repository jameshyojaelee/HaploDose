"""End-to-end recovery and continuous RNA-profile performance candidates."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from haplodose import clones
from haplodose.io import write_json
from haplodose.profile import credible_difference_runs
from haplodose.provenance import start_run, finish_run


def load(name):
    path = Path(__file__).with_name(name+'.py')
    spec = importlib.util.spec_from_file_location('performance_'+name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@contextmanager
def regional_differences(records):
    original = clones.difference_segments
    model = load('regional_amplitude_model')

    def replacement(*args, **kwargs):
        saved = []
        def capture(frame, event, result):
            if frame.f_code is original.__code__ and event == 'return' and result is not None:
                saved.append({key: np.asarray(frame.f_locals[key]).copy() for key in
                    ('d', 'var', 'sub_reset', 'trans', 'init', 'idx', 'post', 'path')})
        assert sys.getprofile() is None
        sys.setprofile(capture)
        try:
            result = original(*args, **kwargs)
        finally:
            sys.setprofile(None)
        assert len(saved) == 1
        payload = saved[0]
        beta, params = args[5:7]
        step = float(np.float32(beta*np.log(1.5)))
        for i, row in enumerate(result):
            if not row.get('testable', True):
                continue
            fitted = model.infer(payload['d'][i], payload['var'][i], payload['sub_reset'],
                payload['trans'], payload['init'], step, step)
            segments = []
            for lo, hi, state in credible_difference_runs(fitted['path'], fitted['post'], payload['sub_reset'], params):
                weights = 1/payload['var'][i, lo:hi+1]
                mean = float(np.sum(weights*payload['d'][i, lo:hi+1])/weights.sum())
                support = fitted['post'][lo:hi+1, state]
                segments.append(dict(first_bin=int(payload['idx'][lo]), last_bin=int(payload['idx'][hi]),
                    n_test_bins=hi-lo+1, direction='higher' if state == 2 else 'lower',
                    mean_diff=mean, z=float(mean*np.sqrt(weights.sum())), posterior=float(support.mean()),
                    min_bin_posterior=float(support.min()), posterior_rule='every_retained_bin'))
            records.append(dict(native_segments=len(row['segments']), regional_segments=len(segments)))
            row['segments'] = segments
        return result

    clones.difference_segments = replacement
    try:
        yield
    finally:
        clones.difference_segments = original


def recovery(args):
    runner = load('membership_stage_panel')
    argv, records = sys.argv, []
    try:
        sys.argv = [str(runner.__file__), '--out', str(args.out/'panel'), '--stage', 'hard',
            '--panel-root', str(args.task.parent/'recovery_20261107/recovery_panel_v1'),
            '--assertion-source', str(args.task.parent/'recovery_20261107/baseline/tests/test_split_mode.py')]
        with regional_differences(records):
            runner.main()
    finally:
        sys.argv = argv
        write_json(records, args.out/'model_calls.json')
    summary = json.loads((args.out/'panel/summary.json').read_text())
    assert len(summary['cases']) == 5 and records
    return dict(cases=summary['cases'], recovered=sum(r['original_assertion_passed'] for r in summary['cases']),
        comparisons=len(records), false_split_power='not yet tested; required before promotion', production_changed=False)


def partition(track, variance, reset):
    output = np.empty(len(track))
    parts = []
    starts = np.flatnonzero(reset)
    assert starts[0] == 0
    for begin, end in zip(starts, np.r_[starts[1:], len(track)]):
        y, w = track[begin:end], 1/variance[begin:end]
        n = len(y)
        penalty = 3*np.log(n) if n > 1 else 0.
        A, B, C = [np.r_[0., np.cumsum(v)] for v in (w, w*y, w*y*y)]
        dp, back = np.full(n+1, np.inf), np.zeros(n+1, int)
        dp[0] = -penalty
        for hi in range(1, n+1):
            cost = (C[hi]-C[:hi])-(B[hi]-B[:hi])**2/(A[hi]-A[:hi])
            values = dp[:hi]+np.maximum(cost, 0)+penalty
            lo = int(np.argmin(values))
            dp[hi], back[hi] = values[lo], lo
        hi = n
        while hi:
            lo = int(back[hi])
            value = float(np.sum(w[lo:hi]*y[lo:hi])/w[lo:hi].sum())
            output[begin+lo:begin+hi] = value
            parts.append(dict(first_bin=int(begin+lo), last_bin=int(begin+hi-1), effect=value))
            hi = lo
    return output, parts


def profile(args):
    root = args.task.parents[1]
    directory = root/f'dev/bulk_amplitude_reference_20261109/refit_{args.refit:03d}_v1'
    base = load('profile')
    source = np.load(directory/'window_scores.npz')
    mapping = pd.read_csv(directory/'gene_map.tsv', sep='\t')
    nwindow = source['dna'].shape[0]
    methods = json.loads((directory/'report.json').read_text())['methods']
    predicted = np.full_like(source['dna'], np.nan)
    records = []
    # All predictions precede DNA scoring or calibration.
    for j, sample in enumerate(source['samples']):
        z = np.load(directory/(str(sample)+'.npz'))
        fitted, segments = partition(z['track'], z['variance'], z['reset'])
        gene_values = np.full(len(mapping), np.nan)
        ok = mapping.shared.to_numpy(bool) & mapping.bin.notna().to_numpy()
        gene_values[ok] = fitted[mapping.loc[ok, 'bin'].to_numpy(int)]
        predicted[:, j] = base.subarm.window_median(gene_values[:, None], mapping.window.to_numpy(), nwindow)[:, 0]
        records.append(dict(sample=str(sample), segments=segments))
    np.save(args.out/'rna_predictions.npy', predicted)
    write_json(records, args.out/'segments.json')
    assert np.isfinite(predicted[source['common']]).all()
    ci = np.flatnonzero(source['split'] == 'calibration')
    rows = []
    for name in [*methods, 'continuous_mdl']:
        values = predicted if name == 'continuous_mdl' else source[name]
        threshold = base.loss_threshold(values, source['neutral'] & source['common'], ci)
        scores = base.score_loss(values, source['loss'], source['neutral'], source['common'], threshold,
            negative=source['gain'] | source['neutral'])
        for j, sample in enumerate(source['samples']):
            ok = source['common'][:, j]
            x, y = values[ok, j], source['dna'][ok, j]
            correlation = float(np.corrcoef(x, y)[0, 1]) if len(x) >= 3 and np.ptp(x) > 0 and np.ptp(y) > 0 else np.nan
            rows.append(dict(sample=str(sample), method=name, split=str(source['split'][j]), pearson=correlation, **scores[j]))
    scores = pd.DataFrame(rows)
    scores.to_csv(args.out/'patient_scores.tsv', sep='\t', index=False)
    scores.groupby(['split', 'method']).agg(pearson=('pearson', 'mean'), recall=('recall_loss', 'mean'),
        fpr=('fpr_loss', 'mean'), precision=('precision_loss', 'mean')).to_csv(args.out/'summary.tsv', sep='\t')
    return dict(methods=[*methods, 'continuous_mdl'], samples=len(source['samples']),
        parameter_search=False, exposed_regression_only=True, production_changed=False)


def ensemble(args):
    base = load('profile')
    collections = {}
    for refit in range(9):
        directory = args.task.parent/f'bulk_amplitude_reference_20261109/refit_{refit:03d}_v1'
        source = np.load(directory/'window_scores.npz')
        names = json.loads((directory/'report.json').read_text())['methods']
        if refit == 0:
            fixed = {k: source[k].copy() for k in ('samples', 'split', 'dna', 'common', 'loss', 'neutral', 'gain')}
        else:
            for key, value in fixed.items():
                assert np.array_equal(source[key], value, equal_nan=True) if value.dtype.kind not in 'US' else np.array_equal(source[key], value)
        for name in names:
            collections.setdefault(name, []).append(source[name].copy())
        collections.setdefault('continuous_mdl', []).append(np.load(args.task/f'profile_sensitivity_v1/refit_{refit:03d}/rna_predictions.npy'))
    predictions = {name: np.median(np.stack(values), axis=0) for name, values in collections.items()}
    np.savez_compressed(args.out/'ensemble_predictions.npz', **predictions)
    ci = np.flatnonzero(fixed['split'] == 'calibration')
    rows = []
    for name, values in predictions.items():
        assert np.isfinite(values[fixed['common']]).all()
        threshold = base.loss_threshold(values, fixed['neutral'] & fixed['common'], ci)
        loss = base.score_loss(values, fixed['loss'], fixed['neutral'], fixed['common'], threshold,
            negative=fixed['gain'] | fixed['neutral'])
        for j, sample in enumerate(fixed['samples']):
            ok = fixed['common'][:, j]
            x, y = values[ok, j], fixed['dna'][ok, j]
            correlation = float(np.corrcoef(x, y)[0, 1]) if len(x) >= 3 and np.ptp(x) > 0 and np.ptp(y) > 0 else np.nan
            rows.append(dict(sample=str(sample), method=name, split=str(fixed['split'][j]), pearson=correlation, **loss[j]))
    scores = pd.DataFrame(rows)
    scores.to_csv(args.out/'patient_scores.tsv', sep='\t', index=False)
    scores.groupby(['split', 'method']).agg(pearson=('pearson', 'mean'), recall=('recall_loss', 'mean'),
        fpr=('fpr_loss', 'mean'), precision=('precision_loss', 'mean')).to_csv(args.out/'summary.tsv', sep='\t')
    return dict(operator='coordinatewise_median_of_all9_reference_predictions',
        all_comparators_ensembled=True, methods=list(predictions), production_changed=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--mode', choices=('recovery', 'profile', 'profile_sensitivity', 'profile_ensemble'), required=True)
    parser.add_argument('--refit', type=int, choices=range(9), default=0)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md')
    if args.mode in ('profile_sensitivity', 'profile_ensemble'):
        inputs['reference_followup'] = args.task/'REFERENCE_FOLLOWUP.md'
    if args.mode == 'profile_ensemble':
        inputs['ensemble_design'] = args.task/'ENSEMBLE_DESIGN.md'
        for refit in range(9):
            inputs[f'candidate{refit}'] = args.task/f'profile_sensitivity_v1/refit_{refit:03d}/rna_predictions.npy'
    if args.mode.startswith('profile'):
        for refit in range(9) if args.mode in ('profile_sensitivity', 'profile_ensemble') else [args.refit]:
            directory = args.task.parent/f'bulk_amplitude_reference_20261109/refit_{refit:03d}_v1'
            for p in sorted(directory.glob('*.npz')):
                inputs[f'profile{refit}:'+p.name] = p
            inputs[f'gene_map{refit}'] = directory/'gene_map.tsv'
            inputs[f'previous_report{refit}'] = directory/'report.json'
    else:
        inputs.update(simulator_manifest=args.task.parent/'recovery_20261107/recovery_panel_v1/manifest.json',
            assertions=args.task.parent/'recovery_20261107/baseline/tests/test_split_mode.py')
    rec = start_run(args.out, inputs, dict(mode=args.mode, refit=args.refit),
        dict(scope='Performance candidate; no promotion before error/power and accuracy gates'), mode='performance_candidate')
    try:
        if args.mode == 'profile_ensemble':
            report = ensemble(args)
        elif args.mode == 'profile_sensitivity':
            reports = []
            for refit in range(9):
                child = args.out/f'refit_{refit:03d}'
                child.mkdir()
                reports.append(profile(argparse.Namespace(task=args.task, out=child, refit=refit)))
            report = dict(reference_refits=9, results=reports, method_promoted=False)
        else:
            report = recovery(args) if args.mode == 'recovery' else profile(args)
        write_json(report, args.out/'report.json')
        finish_run(args.out, rec, backend='CPU full workflow / continuous optimal partition')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
