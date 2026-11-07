"""Content-addressed evidence records; reproducibility is not biological validation."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .io import write_json, guard, _json_default

SCHEMA_VERSION = 2


def file_record(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    before = path.stat()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError(f'Input changed while hashing: {path}')
    return dict(path=str(path), bytes=after.st_size, sha256=digest.hexdigest())


def content_id(record):
    data = json.dumps(record, sort_keys=True, separators=(',', ':'), default=str).encode()
    return hashlib.sha256(data).hexdigest()


def start_run(out, inputs, parameters, contract, *, mode):
    manifest = guard(Path(out) / 'run_manifest.json')
    if manifest.exists():
        raise FileExistsError(f'Refusing to overwrite an existing run manifest: {manifest}')
    package = Path(__file__).resolve().parent
    source = {str(p.relative_to(package)): file_record(p)
              for p in sorted(package.rglob('*')) if p.is_file() and p.suffix in {'.py', '.so', '.tsv'}}
    bundle = dict(files={k: file_record(v) for k, v in inputs.items() if v is not None}, contract=contract)
    bundle['id'] = content_id(bundle)
    execution = dict(input_id=bundle['id'], mode=mode, package_version=__version__, schema_version=2,
                     source=source, parameters=parameters, python=sys.version, executable=sys.executable,
                     package_path=str(package), host=platform.node(), command=sys.argv,
                     environment={d.metadata.get('Name', ''): d.version for d in importlib.metadata.distributions()},
                     allocation={k: os.environ.get(k) for k in ['SLURM_JOB_ID', 'SLURM_JOB_ACCOUNT',
                         'SLURM_JOB_PARTITION', 'SLURM_CPUS_PER_TASK', 'SLURM_MEM_PER_NODE',
                         'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'RAYON_NUM_THREADS']})
    execution['id'] = content_id(execution)
    record = dict(schema_version=2, input=bundle, execution=execution,
                  started=datetime.now(timezone.utc).isoformat(), status='running',
                  evidence=dict(computational_validation='not_assessed_by_run',
                                statistical_identifiability='conditional_relative_RNA',
                                biological_validation='not_established', calibration_status='unvalidated',
                                clone_interpretation='RNA_defined_candidate_groups',
                                tree_interpretation='similarity_dendrogram'))
    # Exclusive creation also protects against two callers starting in the same directory.
    with manifest.open('x') as handle:
        json.dump(record, handle, indent=2, default=_json_default)
        handle.write('\n')
    return record


def finish_run(out, record, *, backend=None, error=None):
    out, record = Path(out), dict(record)
    manifest = out / 'run_manifest.json'
    existing = json.loads(manifest.read_text())
    if (existing.get('execution', {}).get('id') != record.get('execution', {}).get('id')
            or existing.get('started') != record.get('started')
            or existing.get('status') != 'running'):
        raise ValueError('Run manifest does not match this active execution')
    record.update(finished=datetime.now(timezone.utc).isoformat(),
                  status='failed' if error is not None else 'complete', backend=backend)
    if error is not None:
        record['error'] = dict(type=type(error).__name__, message=str(error))
    files = {str(p.relative_to(out)): file_record(p) for p in sorted(out.iterdir())
             if p.is_file() and p.name != 'run_manifest.json' and not p.name.startswith('.')}
    bundle = dict(execution_id=record['execution']['id'], status=record['status'], files=files)
    if 'children' in record:
        bundle['children'] = record['children']
    bundle['id'] = content_id(bundle)
    record['output'] = bundle
    write_json(record, manifest)
    return record
