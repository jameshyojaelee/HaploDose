"""Secondary fixed-seed thinning robustness; execute only after primary recovery.

Simulation draws are technical replicates under one simulator, not independent
biological specimens or evidence of a calibrated specimen false-split rate.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import importlib.metadata
import importlib.util
import inspect
import json
from pathlib import Path
import platform
import shutil
import sys
import time
import traceback

SEEDS = (20261031, 20261107, 20261108, 20261109, 20261110)
ASSERTION = 'test_thin_records_mode_and_recovers_clones'


def record(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def invoke(assertion, out):
    try:
        assertion({'thin': out})
        return dict(passed=True, error=None)
    except Exception:
        return dict(passed=False, error=traceback.format_exc())


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--assertion-source', type=Path, required=True)
    parser.add_argument('--primary-output', type=Path, required=True)
    args = parser.parse_args(argv)
    if not __debug__:
        raise RuntimeError('Python optimization disables original assertions')
    import pandas as pd
    import haplodose
    from haplodose.clones import CloneParams
    from haplodose.pipeline import run
    from haplodose.simulate import write_dataset, simulate, DEFAULT_SPECIMENS

    spec = importlib.util.spec_from_file_location('frozen_split_assertions', args.assertion_source.resolve())
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assertion = getattr(module, ASSERTION)
    source = inspect.getsource(assertion)
    primary = invoke(assertion, args.primary_output.resolve())
    if not primary['passed']:
        raise RuntimeError('Primary original recovery criterion has not passed:\n' + primary['error'])
    args.out.mkdir(parents=True, exist_ok=False)
    package = Path(haplodose.__file__).resolve().parent
    inventory = {d.metadata['Name']: d.version for d in importlib.metadata.distributions() if d.metadata['Name']}
    manifest = dict(evidence_role='secondary_development_simulation_robustness',
                    independent_biological_evaluation=False, simulation_seeds=SEEDS,
                    replication_unit='independent simulation draw, shared simulator and fixed model seed',
                    assertion_source=record(args.assertion_source), assertion_name=ASSERTION,
                    assertion_text_sha256=hashlib.sha256(source.encode()).hexdigest(),
                    primary_gate=primary,
                    primary_files=[record(args.primary_output / name) for name in
                                   ('qc.json', 'cells.tsv.gz', 'differences.tsv')],
                    executed_script=record(__file__), command=sys.argv,
                    environment=dict(python=sys.version, executable=sys.executable, platform=platform.platform(),
                                     packages=inventory, installed_package=str(package)),
                    package_source=[record(p) for p in sorted(package.rglob('*'))
                                    if p.is_file() and '__pycache__' not in p.parts],
                    simulator_defaults=str(inspect.signature(simulate)),
                    default_specimens=DEFAULT_SPECIMENS, params=asdict(CloneParams(split_mode='thin')),
                    chunk=97, cases=[])
    (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    rows = []
    for seed in SEEDS:
        start = time.monotonic()
        case = args.out / f'seed_{seed}'
        case.mkdir()
        entry = dict(seed=seed, passed=False)
        try:
            data = write_dataset(case / 'data', seed=seed)
            original = pd.read_csv(data / 'cells.tsv', sep='\t')
            truth_columns = [c for c in original if c.startswith('truth_')]
            original.drop(columns=truth_columns).to_csv(case / 'model_cells.tsv', sep='\t', index=False)
            entry['input_files'] = [record(p) for p in sorted(data.iterdir()) if p.is_file()]
            entry['model_metadata'] = record(case / 'model_cells.tsv')
            entry['truth_columns_removed_before_fit'] = truth_columns
            run(data / 'counts.h5', case / 'model_cells.tsv', data / 'positions.tsv', case / 'thin',
                chunk=97, params=CloneParams(split_mode='thin'), **module.INPUT_CONTRACT)
            # Truth joins happen only after fitting. Preserve original model outputs.
            fitted = pd.read_csv(case / 'thin/cells.tsv.gz', sep='\t')
            scored = fitted.merge(original[['cell_id'] + truth_columns], on='cell_id',
                                  how='left', validate='one_to_one')
            view = case / 'assertion_view'
            view.mkdir()
            for name in ('qc.json', 'differences.tsv'):
                shutil.copyfile(case / 'thin' / name, view / name)
            scored.to_csv(view / 'cells.tsv.gz', sep='\t', index=False)
            entry.update(invoke(assertion, view))
            qc = json.loads((case / 'thin/qc.json').read_text())
            entry['actual_params'] = qc['params']
            entry['specimens'] = []
            for specimen, frame in scored.groupby('specimen', sort=False):
                if specimen == 'REF':
                    continue
                malignant = frame[frame.malignant_call]
                confusion = pd.crosstab(frame.truth_clone, frame.clone, dropna=False)
                confusion.to_csv(case / f'{specimen}_truth_candidate_confusion.tsv', sep='\t')
                entry['specimens'].append(dict(specimen=specimen, n_cells=len(frame),
                    n_malignant_calls=len(malignant), n_recorded_group_labels=int(frame.clone.nunique()),
                    n_assigned_groups=int(frame.loc[frame.status.eq('assigned'), 'clone'].nunique()),
                    n_malignant_candidate_groups=int(malignant.clone.nunique()),
                    candidate_counts=frame.clone.value_counts(dropna=False).to_dict(),
                    status_counts=frame.status.value_counts(dropna=False).to_dict()))
            # Simulator S3 draws its one biological clone in two ordered depth bands.
            s3 = original[original.specimen.eq('S3')].copy()
            config = next(s for s in DEFAULT_SPECIMENS if s['name'] == 'S3')
            if len(config['clones']) != 1 or config['depth_mix'] != [(0.5, 0.25), (0.5, 1.0)]:
                raise ValueError('Frozen S3 depth-null design changed')
            s3['depth_design_group'] = ['low_0.25x'] * (len(s3) // 2) + ['high_1x'] * (len(s3) - len(s3) // 2)
            s3 = scored.merge(s3[['cell_id', 'depth_design_group']], on='cell_id', validate='one_to_one')
            pd.crosstab(s3.depth_design_group, s3.clone).to_csv(case / 'S3_depth_null_candidate_confusion.tsv', sep='\t')
            entry['s3_depth_null'] = dict(n_cells=len(s3), n_groups=int(s3.clone.nunique()),
                n_malignant_groups=int(s3.loc[s3.malignant_call, 'clone'].nunique()),
                interpretation='One simulated biological clone; depth groups use original draw order, not fitted assignments. '
                               'Descriptive robustness diagnostic, outside the original pass/fail assertion.')
            entry['output_files'] = [record(p) for p in sorted((case / 'thin').iterdir()) if p.is_file()]
        except Exception:
            entry.update(passed=False, execution_error=traceback.format_exc())
        entry['wall_seconds'] = time.monotonic() - start
        (case / 'result.json').write_text(json.dumps(entry, indent=2))
        rows.append(dict(seed=seed, passed=entry['passed'], execution_error=entry.get('execution_error'),
                         assertion_error=entry.get('error'), wall_seconds=entry['wall_seconds']))
        manifest['cases'].append(entry)
        (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=2))
        print(json.dumps(rows[-1]), flush=True)
    pd.DataFrame(rows).to_csv(args.out / 'summary.tsv', sep='\t', index=False)
    manifest['n_passed'] = sum(row['passed'] for row in rows)
    manifest['n_cases'] = len(SEEDS)
    (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
