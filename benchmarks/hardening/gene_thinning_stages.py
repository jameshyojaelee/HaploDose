"""Observational stage trace under the already-frozen known-gene-noise folds."""
import argparse
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
import sys
import traceback

import numpy as np
import pandas as pd
from haplodose import clones, pipeline
from haplodose.io import write_json
from haplodose.provenance import file_record, start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--seed', type=int, choices=(20261108, 20261109, 20261110), required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    arrays, folds = a.out/'arrays', a.out/'folds'
    arrays.mkdir(); folds.mkdir()
    root = Path(__file__).parent
    hook = load('gene_hook', root/'gene_thinning_hook.py')
    scorer = load('original_scorer', root/'recovery_baseline.py')
    src = a.task.parent/f'recovery_20261107/recovery_panel_v1/seed_{a.seed}'
    reference = a.task.parent/f'gene_thinning_20261109/gene_v1/panel/seed_{a.seed}'
    metadata = pd.read_csv(src/'model_cells.tsv', sep='\t')
    assert not any(c.startswith('truth_') for c in metadata)
    s1 = set(np.flatnonzero(metadata.specimen.eq('S1')))
    params = clones.CloneParams(split_mode='thin')
    inputs = dict(script=Path(__file__), hook=root/'gene_thinning_hook.py', scorer=root/'recovery_baseline.py',
        design=a.task/'DESIGN.md', counts=src/'data/counts.h5', metadata=src/'model_cells.tsv',
        positions=src/'data/positions.tsv', previous=reference/'summary.json')
    rec = start_run(a.out, inputs, asdict(params), dict(scope='Observational known-noise development trace, no calibration'),
        mode='gene_thinning_stages')
    targets = {getattr(clones.SpecimenModel, name).__code__: name for name in
        ('_score_candidate', 'event_mixtures', 'agglomerate', 'matched_differences', 'supported', 'merge_unsupported')}
    targets[clones._fit_event_mixture.__code__] = 'event_fit'
    events, errors, active, fold_trace = [], [], {}, []

    def observe(frame, event, result):
        if frame.f_code not in targets or event not in ('call', 'return'):
            return
        try:
            cursor, model = frame, None
            while cursor is not None:
                if isinstance(cursor.f_locals.get('self'), clones.SpecimenModel):
                    model = cursor.f_locals['self']; break
                cursor = cursor.f_back
            if model is None or not set(model.rows.rows).issubset(s1):
                return
            stage, loc = targets[frame.f_code], frame.f_locals
            if event == 'call':
                index = len(events)
                entry = dict(index=index, stage=stage, parent=active.get(id(frame.f_back)),
                    log_start=len(model.log), complete=False)
                events.append(entry); active[id(frame)] = index
                if stage == 'merge_unsupported':
                    path = arrays/f'{index:05d}_merge_input.npz'
                    np.savez_compressed(path, rows=model.rows.rows.copy(), R=loc['R'].copy(),
                        doses=np.asarray(loc['doses']).copy(), pi=loc['pi'].copy())
                    entry['input_arrays'] = str(path.relative_to(a.out))
                return
            index = active.pop(id(frame))
            entry = events[index]
            entry.update(complete=True, log_end=len(model.log))
            saved = dict(rows=model.rows.rows.copy())
            if stage == '_score_candidate':
                saved.update(cells=loc['cells'].copy(), groups=loc['g'].copy())
                entry.update(source=loc['source'], parent_clone=int(loc['k']), accepted=result is not None,
                    members=loc.get('members'), score=None if result is None else float(result[0]))
            elif stage == 'event_fit':
                outer = frame.f_back.f_locals
                if result is None:
                    entry['fitted'] = False
                else:
                    entry.update(fitted=True, gain=float(result['gain']))
                    saved.update(q=result['responsibility'].copy(), parameters=result['parameters'].copy(),
                        cells=outer['cells'].copy(), cols=outer['cols'][outer['part']].copy(),
                        exposure=loc['E'].sum(axis=1))
            elif stage == 'event_mixtures':
                saved['cells'] = loc['cells'].copy()
                entry['returned'] = [dict(source=source, groups=groups.tolist()) for groups, source in result]
            elif stage == 'agglomerate':
                entry.update(context=loc['context'], members=result[0])
            elif stage == 'matched_differences':
                saved['labels'] = loc['labels'].copy()
                entry.update(pairs=loc['pairs'], differences=result)
                parent = frame.f_back
                entry['caller'] = parent.f_code.co_name
                if parent.f_code is clones.SpecimenModel.merge_unsupported.__code__:
                    saved.update(assignment_R=parent.f_locals['assignment_R'].copy(),
                        full_R=parent.f_locals['R'].copy(), doses=np.asarray(parent.f_locals['doses']).copy())
            elif stage == 'supported':
                entry.update(caller=frame.f_back.f_code.co_name, input_segments=loc['segments'], kept=result)
            elif stage == 'merge_unsupported':
                saved['new_doses'] = np.asarray(result[0]).copy()
                entry.update(changed=bool(result[2]), log=model.log[entry['log_start']:entry['log_end']])
            path = arrays/f'{index:05d}_{stage}.npz'
            np.savez_compressed(path, **saved)
            entry['arrays'] = str(path.relative_to(a.out))
        except Exception:
            errors.append(traceback.format_exc())

    try:
        assert sys.getprofile() is None
        with hook.gene_native(fold_trace, folds):
            sys.setprofile(observe)
            try:
                pipeline.run(src/'data/counts.h5', src/'model_cells.tsv', src/'data/positions.tsv', a.out/'fit',
                    params=params, chunk=97, counts_kind='raw', genome_build='hg38',
                    coordinate_system='one-based-closed', gene_namespace='symbol')
            finally:
                sys.setprofile(None)
        # Truth enters only after inference and is not used by any traced operation.
        truth = pd.read_csv(src/'data/cells.tsv', sep='\t')
        fitted = pd.read_csv(a.out/'fit/cells.tsv.gz', sep='\t')
        differences = pd.read_csv(a.out/'fit/differences.tsv', sep='\t')
        pd.testing.assert_frame_equal(fitted, pd.read_csv(reference/'fit/cells.tsv.gz', sep='\t'), check_exact=True)
        pd.testing.assert_frame_equal(differences, pd.read_csv(reference/'fit/differences.tsv', sep='\t'), check_exact=True)
        joined = fitted.merge(truth[['cell_id']+[c for c in truth if c.startswith('truth_')]],
            on='cell_id', validate='one_to_one')
        criteria = scorer.criteria(joined, differences)
        assert criteria == json.loads((reference/'summary.json').read_text())['criteria']
        assert not active and not errors and all(e['complete'] for e in events)
        write_json(events, a.out/'trace.json')
        write_json(fold_trace, a.out/'fold_trace.json')
        write_json(dict(seed=a.seed, criteria=criteria, fit_tables_replay_exactly=True,
            captures=len(events), errors=errors,
            artifacts={str(f.relative_to(a.out)): file_record(f) for f in sorted(arrays.iterdir())}), a.out/'report.json')
    except BaseException as error:
        write_json(dict(events=events, errors=errors), a.out/'partial_trace.json')
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU observational profiling; unchanged known-gene folds')


if __name__ == '__main__':
    main()
