"""Matched original-wheel comparison on the five frozen recovery-panel inputs.

Run with installed-environment Python -I. No simulation generation, selection,
threshold adjustment or primary-pass precondition is performed here.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import platform
import shutil
import sys
import traceback
import zipfile

SEEDS = (20261031, 20261107, 20261108, 20261109, 20261110)


def record(path):
    path = Path(path).resolve()
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=h.hexdigest())


def criteria(cells, differences):
    """Evaluate every original scientific criterion, without short circuiting."""
    import pandas as pd
    s1 = cells[cells.specimen.eq('S1')]
    normal = s1.truth_clone.eq('normal')
    fpr = float(s1.loc[normal, 'malignant_call'].mean())
    tpr = float(s1.loc[~normal, 'malignant_call'].mean())
    s2 = cells[cells.specimen.eq('S2') & cells.malignant_call]
    sub = s1[s1.truth_clone.isin(['clone1', 'clone2']) & s1.malignant_call]
    table = pd.crosstab(sub.truth_clone, sub.clone)
    accuracy = (max(table.values.diagonal().sum(), table.values[::-1].diagonal().sum()) /
                table.values.sum()) if table.values.sum() else None
    d1 = differences[differences.specimen.eq('S1') & ~differences.clone_a.str.endswith('_c0')
                     & ~differences.clone_b.str.endswith('_c0')]
    chromosomes = sorted(set(d1.chromosome))
    return dict(normal_fpr=dict(value=fpr, n=int(normal.sum()), false_calls=int(s1.loc[normal, 'malignant_call'].sum()),
                               criterion='value < 0.1', passed=bool(fpr < .1)),
                tumour_tpr=dict(value=tpr, n=int((~normal).sum()), criterion='value > 0.9', passed=bool(tpr > .9)),
                s2_one_group=dict(value=int(s2.clone.nunique()), criterion='value == 1', passed=bool(s2.clone.nunique() == 1)),
                s1_two_groups=dict(value=int(sub.clone.nunique()), criterion='value == 2', passed=bool(sub.clone.nunique() == 2)),
                assignment_accuracy=dict(value=None if accuracy is None else float(accuracy),
                    n=int(table.values.sum()), criterion='value > 0.65',
                    passed=bool(accuracy is not None and accuracy > .65),
                    interpretation='Original best two-row diagonal accuracy; interpret with separate two-group criterion'),
                difference_chromosomes=dict(value=chromosomes, criterion="set(value) == {'chr12'}", passed=chromosomes == ['chr12']))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--panel-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--assertion-source', type=Path, required=True)
    parser.add_argument('--wheel', type=Path, required=True)
    args = parser.parse_args(argv)
    if not sys.flags.isolated or not __debug__:
        raise RuntimeError('Use Python -I without -O: isolated installed imports and active assertions are required')
    import pandas as pd
    import haplodose
    from haplodose.clones import CloneParams
    from haplodose.pipeline import run
    package = Path(haplodose.__file__).resolve().parent
    wheel_binding = []
    with zipfile.ZipFile(args.wheel) as archive:
        members = [name for name in archive.namelist() if name.startswith('haplodose/') and not name.endswith('/')]
        if not members:
            raise ValueError('Wheel has no haplodose package members')
        for name in members:
            installed = package / Path(name).relative_to('haplodose')
            expected = hashlib.sha256(archive.read(name)).hexdigest()
            observed = record(installed)
            if observed['sha256'] != expected:
                raise ValueError(f'Installed package does not match required baseline wheel: {name}')
            wheel_binding.append(observed)
    spec = importlib.util.spec_from_file_location('baseline_assertions', args.assertion_source.resolve())
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assertion = module.test_thin_records_mode_and_recovers_clones
    source_manifest = json.loads((args.panel_root / 'manifest.json').read_text())
    if tuple(source_manifest['simulation_seeds']) != SEEDS:
        raise ValueError('Original five-seed panel manifest changed')
    args.out.mkdir(parents=True, exist_ok=False)
    manifest = dict(role='matched_baseline_secondary_development_simulation', seeds=SEEDS,
        package=str(package), script=record(__file__), assertion_source=record(args.assertion_source),
        wheel=record(args.wheel), installed_wheel_member_binding=wheel_binding,
        original_panel_manifest=record(args.panel_root / 'manifest.json'),
        params=asdict(CloneParams(split_mode='thin')), chunk=97,
        environment=dict(python=sys.version, executable=sys.executable, isolated=bool(sys.flags.isolated),
            platform=platform.platform(), packages={d.metadata['Name']: d.version
                for d in importlib.metadata.distributions() if d.metadata['Name']}),
        command=sys.argv, interpretation='Same simulation inputs and thresholds; five development simulation draws, no biological replication', cases=[])
    rows = []
    for seed in SEEDS:
        source = args.panel_root / f'seed_{seed}'
        case = args.out / f'seed_{seed}'
        case.mkdir()
        entry = dict(seed=seed)
        try:
            previous = next(r for r in source_manifest['cases'] if r['seed'] == seed)
            paths = [source / 'data/counts.h5', source / 'data/cells.tsv', source / 'data/positions.tsv', source / 'model_cells.tsv']
            entry['inputs'] = [record(p) for p in paths]
            old_records = previous['input_files'] + [previous['model_metadata']]
            expected_by_path = {str(Path(r['path']).resolve()): r['sha256'] for r in old_records}
            for r in entry['inputs']:
                if expected_by_path.get(r['path']) != r['sha256']:
                    raise ValueError('Saved panel input no longer matches original panel hash: ' + r['path'])
            metadata = pd.read_csv(source / 'model_cells.tsv', sep='\t')
            if any(c.startswith('truth_') for c in metadata):
                raise ValueError('Truth must be excluded from fitting metadata')
            run(source / 'data/counts.h5', source / 'model_cells.tsv', source / 'data/positions.tsv',
                case / 'thin', chunk=97, params=CloneParams(split_mode='thin'), **module.INPUT_CONTRACT)
            fitted = pd.read_csv(case / 'thin/cells.tsv.gz', sep='\t')
            original = pd.read_csv(source / 'data/cells.tsv', sep='\t')
            truth = [c for c in original if c.startswith('truth_')]
            scored = fitted.merge(original[['cell_id'] + truth], on='cell_id', how='left', validate='one_to_one')
            view = case / 'assertion_view'
            view.mkdir()
            for name in ('qc.json', 'differences.tsv'):
                shutil.copyfile(case / 'thin' / name, view / name)
            scored.to_csv(view / 'cells.tsv.gz', sep='\t', index=False)
            try:
                assertion({'thin': view})
                entry['original_assertion_passed'] = True
            except Exception:
                entry.update(original_assertion_passed=False, assertion_error=traceback.format_exc())
            differences = pd.read_csv(case / 'thin/differences.tsv', sep='\t')
            entry['criteria'] = criteria(scored, differences)
            candidate_cells = pd.read_csv(source / 'assertion_view/cells.tsv.gz', sep='\t')
            candidate_differences = pd.read_csv(source / 'thin/differences.tsv', sep='\t')
            entry['candidate_criteria'] = criteria(candidate_cells, candidate_differences)
            entry['candidate_original_assertion_passed'] = previous['passed']
            entry['candidate_scoring_sources'] = [record(source / name) for name in
                ('assertion_view/cells.tsv.gz', 'thin/differences.tsv', 'result.json')]
            entry['split_mode_recorded'] = json.loads((case / 'thin/qc.json').read_text())['params']['split_mode'] == 'thin'
            for specimen, frame in scored.groupby('specimen', sort=False):
                pd.crosstab(frame.truth_clone, frame.clone, dropna=False).to_csv(case / f'{specimen}_confusion.tsv', sep='\t')
            entry['outputs'] = [record(p) for p in sorted((case / 'thin').iterdir()) if p.is_file()]
        except Exception:
            entry.update(original_assertion_passed=False, execution_error=traceback.format_exc())
        (case / 'result.json').write_text(json.dumps(entry, indent=2))
        row = dict(seed=seed, original_assertion_passed=entry['original_assertion_passed'],
                   execution_error=entry.get('execution_error'))
        for name, result in entry.get('criteria', {}).items():
            row[name + '_passed'] = result['passed']
            row[name + '_value'] = result['value']
            comparison = entry.get('candidate_criteria', {}).get(name)
            if comparison is None:
                continue
            row['candidate_' + name + '_passed'] = comparison['passed']
            row['candidate_' + name + '_value'] = comparison['value']
            row[name + '_introduced_failure'] = bool(result['passed'] and not comparison['passed'])
            row[name + '_recovered_failure'] = bool(not result['passed'] and comparison['passed'])
        rows.append(row)
        manifest['cases'].append(entry)
        (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=2))
        print(json.dumps(row), flush=True)
    pd.DataFrame(rows).to_csv(args.out / 'summary.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
