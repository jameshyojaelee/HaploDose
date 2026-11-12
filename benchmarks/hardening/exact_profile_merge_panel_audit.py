"""Audit exposed full-workflow regression without reinterpreting it as validation."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    path = Path(__file__).with_name('gene_veto_fresh_audit.py')
    spec = importlib.util.spec_from_file_location('independent_panel_audit', path)
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    paths = sorted(a.task.glob('panel_*_v1/results.json'))
    assert len(paths)==12
    rec = start_run(a.out, dict(script=Path(__file__), independent_scorer=path,
        design=a.task/'PANEL_DESIGN.md', **{str(i):f for i,f in enumerate(paths)}), {},
        dict(scope='Exposed-panel regression only; conditional bounds descriptive, biological n=0'), mode='exact_merge_panel_audit')
    try:
        rows = [r for path in paths for r in json.loads(path.read_text())]
        assert len(rows)==600
        keys = [(r['condition'], r['rep'], r['units']) for r in rows]
        assert len(set(keys))==600
        assert set(keys)=={(c,i,u) for c in ('matched','state','two_clone') for i in range(100) for u in ('raw','effective')}
        checked, inputs, folds = [], {}, 0
        for row in rows:
            old = Path(row['prior']['case'])
            folder = Path(row['out'])
            truth = pd.read_csv(old/'input/truth_cells.tsv', sep='\t')
            expected = audit.independent_score(folder/'fit', truth)
            for name, value in expected.items():
                assert row[name]==value, (row['condition'], row['rep'], row['units'], name)
            before = audit.independent_score(old/'candidate/fit', truth)
            for name, value in before.items():
                assert row['prior'][name]==value
            original = json.loads((old/'candidate/fit/run_manifest.json').read_text())['input']
            current = json.loads((folder/'fit/run_manifest.json').read_text())['input']
            assert current==original
            for record in current['files'].values():
                path = record['path']
                if path not in inputs:
                    inputs[path] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
                assert inputs[path]==record['sha256']
            for fold in (folder/'folds').glob('*.npz'):
                with np.load(fold) as x, np.load(old/'candidate/folds'/fold.name) as y:
                    assert set(x.files)==set(y.files)
                    assert all(np.array_equal(x[k],y[k]) for k in x.files)
                folds += 1
            checked.append(dict(condition=row['condition'], rep=row['rep'], units=row['units'],
                                before=before, after=expected))
        summaries = []
        for condition in ('matched','state','two_clone'):
            for units in ('raw','effective'):
                selected = [r for r in checked if r['condition']==condition and r['units']==units]
                upper = condition!='two_clone'
                event = 'null_error' if upper else 'exact_chr12_recovery'
                k = sum(r['after'][event] for r in selected)
                b, independent = audit.bound(k,100,upper,.05/6)
                summaries.append(dict(condition=condition,units=units,events=k,n=100,
                    bound=b,bound_direction='upper' if upper else 'lower',
                    independent_binomial_inversion=independent,
                    old_events=sum(r['before'][event] for r in selected),
                    mean_all_cell_assignment=sum(r['after']['assignment_correct'] for r in selected)/60000))
        # Retain the entire historical experiment, including its comparator.
        frozen_files = {}
        previous = a.task.parent/'gene_veto_fresh_20261109'
        for report in sorted(previous.glob('chunk_*_v1/results.json')):
            for row in json.loads(report.read_text()):
                record = row['fit_manifest']
                path = Path(record['path'])
                assert hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256']
                manifest = json.loads(path.read_text())
                for output in manifest['output']['files'].values():
                    filename = output['path']
                    if filename not in frozen_files:
                        h = hashlib.sha256()
                        with Path(filename).open('rb') as stream:
                            for block in iter(lambda:stream.read(1048576), b''):
                                h.update(block)
                        frozen_files[filename] = h.hexdigest()
                    assert frozen_files[filename]==output['sha256'], filename
        write_json(dict(all_checks_passed=True, fits_checked=600, input_records_checked=1800,
            unique_files_hashed=len(inputs), exact_folds_checked=folds, summaries=summaries,
            frozen_original_output_files_verified=len(frozen_files), frozen_original_fit_manifests_verified=600,
            cases=checked, biological_specimens=0, scientific_criteria_changed=False,
            interpretation='Retrospective exposed-panel regression; no prospective or biological guarantee'), a.out/'report.json')
        finish_run(a.out, rec, backend='CPU exhaustive public-output matching; exact inputs/folds; beta vs binomial inversion')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
