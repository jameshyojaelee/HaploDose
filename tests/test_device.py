"""Device selection (haplodose.device): CPU default, GPU opt in, clear errors. Each case runs in a fresh process,
because JAX binds to a platform once per process."""
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

PY = [sys.executable, '-s', '-c']


def fresh(code, **env):
    e = dict(os.environ, **env)
    e.pop('JAX_PLATFORMS', None)
    return subprocess.run(PY + [code], capture_output=True, text=True, env=e, timeout=600)


def gpu_available():
    r = fresh('import jax; print(any(d.platform == "gpu" for d in jax.devices()))')
    return r.returncode == 0 and r.stdout.strip().endswith('True')


def test_import_does_not_start_jax_and_cpu_is_default():
    r = fresh('import haplodose.cli, haplodose.pipeline, haplodose.bulk\n'
              'from jax._src import xla_bridge\n'
              'started = xla_bridge.backends_are_initialized()\n'
              'from haplodose.device import select\n'
              'print(started, select())')
    assert r.returncode == 0, r.stderr
    assert r.stdout.split()[-2:] == ['False', 'cpu']


def test_invalid_device_rejected():
    from haplodose.device import select
    with pytest.raises(ValueError):
        select('tpu')


def test_gpu_without_visible_gpu_fails_clearly():
    r = fresh('from haplodose.device import select; select("gpu")', CUDA_VISIBLE_DEVICES='')
    assert r.returncode != 0
    assert 'device gpu' in r.stderr


def test_auto_falls_back_to_cpu_without_gpu():
    r = fresh('from haplodose.device import select; print(select("auto"))', CUDA_VISIBLE_DEVICES='')
    assert r.returncode == 0, r.stderr
    assert r.stdout.split()[-1] == 'cpu'


def test_device_cannot_change_after_jax_started():
    r = fresh('from haplodose.device import select; select("cpu"); select("gpu")')
    assert r.returncode != 0
    assert 'already running on cpu' in r.stderr


RUN = '''
import json, sys
from pathlib import Path
from haplodose.pipeline import run
data, out, device = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
qc = run(data / 'counts.h5', data / 'cells.tsv', data / 'positions.tsv', out, chunk=97, device=device,
         counts_kind='raw', genome_build='hg38', coordinate_system='one-based-closed', gene_namespace='symbol')
print(json.dumps(dict(device=qc['device'], backend=qc['backend'])))
'''


@pytest.mark.skipif(not gpu_available(), reason='no GPU visible')
def test_gpu_matches_cpu(tmp_path):
    """Same simulated data on CPU and GPU: the clone calls and profiles must agree closely. Float reduction order
    and the GPU HMM algorithm (associative scan) differ, so exact equality is not required."""
    from sklearn.metrics import adjusted_rand_score
    from haplodose.simulate import write_dataset
    data = write_dataset(tmp_path / 'data')
    outs = {}
    for device in ('cpu', 'gpu'):
        r = subprocess.run([sys.executable, '-s', '-c', RUN, str(data), str(tmp_path / device), device],
                           capture_output=True, text=True, timeout=3600)
        assert r.returncode == 0, r.stderr[-2000:]
        assert json.loads(r.stdout.strip().splitlines()[-1]) == dict(device=device, backend=device)
        outs[device] = tmp_path / device
    cells = {d: pd.read_csv(o / 'cells.tsv.gz', sep='\t').set_index('cell_id') for d, o in outs.items()}
    a, b = cells['cpu'], cells['gpu'].reindex(cells['cpu'].index)
    assert (a.malignant_call == b.malignant_call).mean() > 0.99
    for spec, c in a.groupby('specimen'):
        assert adjusted_rand_score(c.clone.astype(str), b.loc[c.index, 'clone'].astype(str)) > 0.95, spec
    prof = {d: pd.read_csv(o / 'profiles.tsv.gz', sep='\t') for d, o in outs.items()}
    m = prof['cpu'].merge(prof['gpu'], on=['specimen', 'clone', 'bin'], suffixes=('_cpu', '_gpu'))
    assert len(m) > 0.95 * len(prof['cpu'])
    assert np.corrcoef(m.log_ratio_cpu, m.log_ratio_gpu)[0, 1] > 0.999
