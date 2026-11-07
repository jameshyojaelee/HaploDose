"""Observe the unchanged default thin-recovery fixture on a compute node.

Python return-event profiling copies existing intermediates; no numerical operation,
random stream, parameter, candidate or acceptance rule is replaced.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
import sys
import traceback

import numpy as np
import pandas as pd

from haplodose import clones, pipeline, profile
from haplodose.io import new_output_dir, write_json
from haplodose.provenance import file_record, start_run, finish_run
from haplodose.simulate import write_dataset


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--assertion-source', type=Path, required=True,
                    help='Frozen unchanged tests/test_split_mode.py')
    ap.add_argument('--truth-diagnostics', action='store_true')
    args = ap.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    out = new_output_dir(args.out)
    data = write_dataset(out / 'data')  # Original default simulation, no kwargs.
    meta = pd.read_csv(data / 'cells.tsv', sep='\t')
    s1rows = set(np.flatnonzero(meta.specimen.eq('S1')))
    params = clones.CloneParams(split_mode='thin')
    inputs = {p.name: p for p in sorted(data.iterdir()) if p.is_file()}
    inputs.update(runner=Path(__file__), original_assertions=args.assertion_source.resolve())
    record = start_run(out, inputs, asdict(params),
        dict(counts_kind='raw', genome_build='hg38', coordinate_system='one-based-closed',
             gene_namespace='symbol', diagnostic_scope='unchanged default simulation S1; no truth used in fit'),
        mode='thin_recovery_trace')
    arrays = out / 'arrays'
    arrays.mkdir()
    events, errors = [], []
    targets = {clones.SpecimenModel._score_candidate.__code__: 'candidate',
               clones.SpecimenModel.differences.__code__: 'differences',
               clones.SpecimenModel.event_mixtures.__code__: 'event_mixtures',
               profile.difference_segments.__code__: 'difference_segments'}

    def observe(frame, event, result):
        if event != 'return' or frame.f_code not in targets:
            return
        try:
            cursor, model, context = frame, None, {}
            while cursor is not None:
                loc = cursor.f_locals
                if model is None and isinstance(loc.get('self'), clones.SpecimenModel):
                    model = loc['self']
                for key in ['source', 'context', 'k']:
                    if key in loc and key not in context and isinstance(loc[key], (str, int, np.integer)):
                        context[key] = int(loc[key]) if isinstance(loc[key], np.integer) else loc[key]
                cursor = cursor.f_back
            if model is None or not len(model.rows.rows) or not set(model.rows.rows).issubset(s1rows):
                return
            label = targets[frame.f_code]
            entry = dict(index=len(events), stage=label, context=context)
            loc = frame.f_locals
            keys = {'candidate': ['cells', 'g', 'Yg', 'Lg', 'XLg'],
                    'event_mixtures': ['Yc', 'Lc', 'xc', 'lr', 'uw', 'cells', 'cols', 'Yp',
                                       'Lp', 'Ep', 'Rnb', 'r_bin', 'chrom', 's'],
                    'differences': ['Ya', 'La', 'XLa', 'Yb', 'Lb', 'XLb', 'shift'],
                    'difference_segments': ['Ya', 'La', 'Yb', 'Lb', 'reset', 'idx', 'sub_reset',
                                            'shift', 'd', 'v', 'tau', 'var', 'post', 'path', 'testable']}[label]
            saved = {key: np.asarray(loc[key]).copy() for key in keys if key in loc and loc[key] is not None}
            saved['global_cell_rows'] = model.rows.rows.copy()
            if label in ['candidate', 'event_mixtures']:
                entry['clone_parameters'] = asdict(model.params)
                for key in ['phi', 'thin_size', 'gamma', 'base', 'reset', 'propose_mask', 'test_mask']:
                    if hasattr(model, key):
                        saved['model_' + key] = np.asarray(getattr(model, key)).copy()
            if label == 'candidate':
                entry.update(accepted=result is not None, group_sizes=np.bincount(loc['g']).tolist())
                if result is not None:
                    entry['score'] = float(result[0])
            elif label == 'event_mixtures':
                # Preserve every ranked region before deduplication and all returned splits.
                # These are existing optimizer outputs, not independently recomputed fits.
                entry['found'] = []
                for i, (gain, groups, source) in enumerate(loc.get('found', [])):
                    saved[f'found_groups_{i}'] = np.asarray(groups).copy()
                    entry['found'].append(dict(gain=float(gain), source=source,
                                              group_sizes=np.bincount(groups).tolist()))
                entry['returned'] = []
                for i, (groups, source) in enumerate(result):
                    saved[f'returned_groups_{i}'] = np.asarray(groups).copy()
                    entry['returned'].append(dict(source=source, group_sizes=np.bincount(groups).tolist()))
            elif label == 'difference_segments':
                entry.update(results=result, beta=float(loc['beta']), parameters=asdict(loc['params']))
            else:
                entry['results'] = result
            # All S1 comparisons are small in the original fixed fixture; preserve every pair.
            path = arrays / f'{len(events):04d}_{label}.npz'
            np.savez_compressed(path, **saved)
            entry['arrays'] = str(path.relative_to(out))
            events.append(entry)
        except Exception:
            errors.append(traceback.format_exc())

    try:
        if sys.getprofile() is not None:
            raise RuntimeError('An existing profiler is active; refusing to replace it')
        sys.setprofile(observe)
        try:
            pipeline.run(data / 'counts.h5', data / 'cells.tsv', data / 'positions.tsv', out / 'fit',
                         chunk=97, params=params, counts_kind='raw', genome_build='hg38',
                         coordinate_system='one-based-closed', gene_namespace='symbol')
        finally:
            sys.setprofile(None)
        spec = importlib.util.spec_from_file_location('_preserved_thin_assertions', args.assertion_source)
        assertions = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(assertions)
        gate = dict(status='passed', assertion='test_thin_records_mode_and_recovers_clones')
        try:
            assertions.test_thin_records_mode_and_recovers_clones({'thin': out / 'fit'})
        except AssertionError:
            gate.update(status='failed', traceback=traceback.format_exc())
        write_json(events, out / 'trace.json')
        write_json(dict(gate=gate, capture_errors=errors, captures=len(events),
            instrumentation='Return-event local-variable copies only; no numerical recomputation',
            arrays={str(p.relative_to(out)): file_record(p) for p in sorted(arrays.iterdir())}), out / 'summary.json')
        if args.truth_diagnostics:
            # Truth is first accessed after fit and unchanged assertions, never to select candidates.
            cells = pd.read_csv(out / 'fit/cells.tsv.gz', sep='\t')
            sub = cells[cells.specimen.eq('S1')]
            pd.crosstab(sub.truth_clone, sub.clone, dropna=False).to_csv(out / 'truth_crosstab.tsv', sep='\t')
        if errors:
            raise RuntimeError('Diagnostic captures failed; see summary.json')
        finish_run(out, record, backend='cpu')
        print(json.dumps(gate))
    except Exception as error:
        finish_run(out, record, backend='cpu', error=error)
        raise


if __name__ == '__main__':
    main()
