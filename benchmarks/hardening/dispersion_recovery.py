"""Paired all-five original recovery cases, changing only dispersion latent weights."""
import argparse
from contextlib import nullcontext
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
import shutil
import time

import jax
import pandas as pd
from haplodose.clones import CloneParams
from haplodose.pipeline import run
from haplodose.provenance import start_run, finish_run, file_record


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--arm', choices=['MAP', 'soft'], required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    if jax.default_backend() != 'cpu':
        raise RuntimeError('CPU backend required inside inference')
    panel = a.task.parent/'recovery_20261107/recovery_panel_v1'
    assertion_path = a.task.parent/'recovery_20261107/baseline/tests/test_split_mode.py'
    assertions = load('assertions', assertion_path)
    score = load('score', Path(__file__).with_name('recovery_baseline.py'))
    hook = load('dispersion', Path(__file__).with_name('dispersion_weighted.py'))
    params = CloneParams(split_mode='thin')
    rec = start_run(a.out, dict(script=Path(__file__), design=a.task/'RECOVERY_DESIGN.md',
        hook=Path(hook.__file__), scorer=Path(score.__file__), assertions=assertion_path,
        panel_manifest=panel/'manifest.json', exact_audit=a.task/'matched_audit_v1/run_manifest.json'),
        dict(arm=a.arm, params=asdict(params)), dict(scope='Original exposed recovery panel, no promotion'),
        mode='dispersion_recovery')
    try:
        audit = json.loads((a.task/'matched_audit_v1/run_manifest.json').read_text())
        assert audit['status'] == 'complete'
        manifest = json.loads((panel/'manifest.json').read_text())
        results = []
        for seed in score.SEEDS:
            src = panel/f'seed_{seed}'
            folder = a.out/f'seed_{seed}'
            folder.mkdir()
            start = time.monotonic()
            entry = next(row for row in manifest['cases'] if row['seed'] == seed)
            for record in entry['input_files']+[entry['model_metadata']]:
                assert file_record(record['path'])['sha256'] == record['sha256']
            metadata = pd.read_csv(src/'model_cells.tsv', sep='\t')
            assert not any(c.startswith('truth_') for c in metadata)
            trace = []
            with (hook.weighted_dispersion(trace) if a.arm == 'soft' else nullcontext()):
                run(src/'data/counts.h5', src/'model_cells.tsv', src/'data/positions.tsv', folder/'fit',
                    chunk=97, params=params, **assertions.INPUT_CONTRACT)
            truth = pd.read_csv(src/'data/cells.tsv', sep='\t')
            fitted = pd.read_csv(folder/'fit/cells.tsv.gz', sep='\t')
            truth_columns = [c for c in truth if c.startswith('truth_')]
            joined = fitted.merge(truth[['cell_id']+truth_columns], on='cell_id', validate='one_to_one')
            view = folder/'assertion_view'
            view.mkdir()
            joined.to_csv(view/'cells.tsv.gz', sep='\t', index=False)
            for name in ('qc.json', 'differences.tsv'):
                shutil.copyfile(folder/'fit'/name, view/name)
            passed = True
            try:
                assertions.test_thin_records_mode_and_recovers_clones({'thin': view})
            except AssertionError:
                passed = False
            checks = score.criteria(joined, pd.read_csv(folder/'fit/differences.tsv', sep='\t'))
            row = dict(seed=seed, arm=a.arm, original_assertion_passed=passed, criteria=checks,
                seconds=time.monotonic()-start, s3_groups=int(joined.loc[joined.specimen.eq('S3') & joined.malignant_call, 'clone'].nunique()))
            results.append(row)
            (folder/'summary.json').write_text(json.dumps(row, indent=2)+'\n')
            (folder/'dispersion_trace.json').write_text(json.dumps(trace, indent=2)+'\n')
            print(json.dumps(row), flush=True)
        (a.out/'summary.json').write_text(json.dumps(dict(arm=a.arm, cases=results, biological_specimens=0,
            production_changed=False, scientific_criteria_changed=False), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU verified inside full raw-count inference')


if __name__ == '__main__':
    main()
