"""Replay locked merge mechanisms; retain all original public-output criteria."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose import clones, pipeline
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@contextmanager
def merge_calls(records, directory):
    """Profile calls without replacing the method used by stage-hook identity checks."""
    assert sys.getprofile() is None
    pending = {}
    code = clones.SpecimenModel.merge_unsupported.__code__
    def observe(frame, event, value):
        if frame.f_code is not code:
            return
        if event == 'call':
            loc = frame.f_locals
            pending[id(frame)] = dict(R_before=loc['R'].copy(), pi_before=loc['pi'].copy(),
                                      doses_before=np.asarray(loc['doses']).copy())
        elif event == 'return' and value is not None:
            arrays = pending.pop(id(frame))
            doses, pi, merged = value
            arrays.update(doses_after=np.asarray(doses), pi_after=pi, R_after=frame.f_locals['R'])
            path = directory/f'merge_{len(records):03d}.npz'
            np.savez_compressed(path, **arrays)
            records.append(dict(arrays=str(path), merged=bool(merged),
                                statistical_merge=bool(frame.f_locals['statistical_merge'])))
    sys.setprofile(observe)
    try:
        yield
    finally:
        sys.setprofile(None)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    trace = load('merge_trace', root/'gene_state_trace.py')
    gene = load('merge_gene', root/'gene_thinning_hook.py')
    stage = load('merge_stage', root/'membership_weighted_stages.py')
    veto = load('merge_veto', root/'gene_veto_weights.py')
    scoring = load('merge_scoring', root/'gene_veto_fresh.py')
    prior = a.task.parent/'gene_state_trace_20261109'
    cases = json.loads((prior/'trace_v1/report.json').read_text())
    rec = start_run(a.out, dict(script=Path(__file__), design=a.task/'DESIGN.md',
        prior_audit=prior/'audit_v1/report.json'), {},
        dict(scope='Six exposed regression mechanisms, biological n=0'), mode='exact_profile_merge')
    summaries = []
    try:
        assert json.loads((prior/'audit_v1/report.json').read_text())['all_checks_passed']
        for selected in cases:
            seed, condition = selected['seed'], selected['condition']
            case = Path(selected['prior'])
            out = a.out/condition/f'seed_{seed}'
            out.mkdir(parents=True)
            for name in ('arrays', 'folds', 'veto'):
                (out/name).mkdir()
            stages, gt, st, vt, merges = [], [], [], [], []
            with gene.gene_native(gt, out/'folds', case_seeds=(seed,)), \
                    stage.weighted_stages('both', st, record_memberships=True), \
                    veto.observe_weights(vt, out/'veto', apply=True, weight_units='raw'), \
                    trace.capture(stages, out/'arrays'), merge_calls(merges, out/'arrays'):
                pipeline.run(case/'input/counts.h5', case/'input/model_cells.tsv',
                    case/'input/positions.tsv', out/'fit', params=clones.CloneParams(split_mode='thin'),
                    chunk=97, counts_kind='raw', genome_build='hg38',
                    coordinate_system='one-based-closed', gene_namespace='symbol')
            write_json(dict(stages=stages, gene=gt, weighted=st, veto=vt, merges=merges), out/'trace.json')
            uninstrumented = a.task/'mechanisms_v1'/condition/f'seed_{seed}'/'fit'
            if uninstrumented.exists() and uninstrumented != out/'fit':
                for name in ('cells.tsv.gz', 'clones.tsv', 'differences.tsv', 'segments.tsv', 'profiles.tsv.gz'):
                    pd.testing.assert_frame_equal(pd.read_csv(out/'fit'/name, sep='\t'),
                        pd.read_csv(uninstrumented/name, sep='\t'), check_exact=True)
                assert np.array_equal(np.load(out/'fit/posteriors.npy'),
                    np.load(uninstrumented/'posteriors.npy'), equal_nan=True)
            metadata = pd.read_csv(case/'input/truth_cells.tsv', sep='\t')
            before = scoring.score(case/'candidate/fit', metadata, condition)
            after = scoring.score(out/'fit', metadata, condition)
            # Same count-splitting realizations isolate the production change.
            for fold in (out/'folds').glob('*.npz'):
                old = Path(selected['out'])/'folds'/fold.name
                x, y = np.load(fold), np.load(old)
                assert set(x.files) == set(y.files)
                assert all(np.array_equal(x[k], y[k]) for k in x.files)
            identical = {}
            for name in ('cells.tsv.gz', 'clones.tsv', 'differences.tsv', 'segments.tsv', 'profiles.tsv.gz'):
                x, y = (pd.read_csv(folder/name, sep='\t') for folder in (case/'candidate/fit', out/'fit'))
                identical[name] = bool(x.equals(y))
            summaries.append(dict(condition=condition, rep=selected['rep'], seed=seed,
                before=before, after=after, unchanged_tables=identical,
                out=str(out), trace=file_record(out/'trace.json')))
            write_json(summaries, a.out/'report.json')
        finish_run(a.out, rec, backend='CPU; new isolated wheel; exposed mechanism replay')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
