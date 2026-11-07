"""Independent hashes and failure states for the schema-2 execution chain."""
import hashlib
import gzip
import json

import h5py
import numpy as np
import pandas as pd
import pytest

from haplodose import provenance


CONTRACT = dict(counts_kind='raw', genome_build='hg38', coordinate_system='one-based-closed',
                gene_namespace='symbol')


def test_compressed_tables_have_no_wall_clock_dependency(tmp_path):
    from haplodose.io import write_table
    paths = []
    for name in ('first', 'second'):
        directory = tmp_path / name
        directory.mkdir()
        path = directory / 'profile.tsv.gz'
        write_table(pd.DataFrame({'gene': ['A'], 'count': [17]}), path)
        with gzip.open(path, 'rt') as handle:
            assert handle.read() == 'gene\tcount\nA\t17\n'
        assert int.from_bytes(path.read_bytes()[4:8], 'little') == 0
        paths.append(path)
    assert paths[0].read_bytes() == paths[1].read_bytes()


def expected_id(record):
    payload = {key: value for key, value in record.items() if key != 'id'}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()


def test_manifest_links_exact_input_execution_and_output_content(tmp_path):
    data = tmp_path / 'input.tsv'
    data.write_bytes(b'gene\tcount\nG1\t4\n')
    out = tmp_path / 'out'
    out.mkdir()
    initial = provenance.start_run(out, {'counts': data, 'absent_optional': None}, {'seed': 11}, CONTRACT,
                                   mode='test')
    assert initial['input']['files']['counts']['sha256'] == hashlib.sha256(data.read_bytes()).hexdigest()
    assert initial['input']['id'] == expected_id(initial['input'])
    assert initial['execution']['id'] == expected_id(initial['execution'])
    assert initial['execution']['input_id'] == initial['input']['id']
    (out / 'profiles.tsv').write_bytes(b'bin\tvalue\n0\t0.2\n')
    finished = provenance.finish_run(out, initial, backend='cpu')
    assert finished['status'] == 'complete'
    assert finished['output']['execution_id'] == initial['execution']['id']
    assert finished['output']['id'] == expected_id(finished['output'])
    assert finished['output']['files']['profiles.tsv']['sha256'] == hashlib.sha256(
        (out / 'profiles.tsv').read_bytes()).hexdigest()
    assert 'run_manifest.json' not in finished['output']['files']
    assert finished['evidence']['biological_validation'] == 'not_established'


def test_restarting_or_stale_finishing_cannot_overwrite_evidence(tmp_path):
    data = tmp_path / 'counts'
    data.write_text('1\n')
    out = tmp_path / 'out'
    out.mkdir()
    record = provenance.start_run(out, {'counts': data}, {}, CONTRACT, mode='test')
    manifest = out / 'run_manifest.json'
    original = manifest.read_bytes()
    with pytest.raises(FileExistsError):
        provenance.start_run(out, {'counts': data}, {}, CONTRACT, mode='test')
    assert manifest.read_bytes() == original
    stale = dict(record, execution=dict(record['execution'], id='different-execution'))
    with pytest.raises(ValueError, match='active execution'):
        provenance.finish_run(out, stale, backend='cpu')
    assert manifest.read_bytes() == original
    provenance.finish_run(out, record, backend='cpu')
    completed = manifest.read_bytes()
    with pytest.raises(ValueError, match='active execution'):
        provenance.finish_run(out, record, error=RuntimeError('late failure'))
    assert manifest.read_bytes() == completed


def test_failed_run_records_partial_artifacts_as_failed(tmp_path):
    data = tmp_path / 'counts'
    data.write_text('1\n')
    out = tmp_path / 'out'
    out.mkdir()
    record = provenance.start_run(out, {'counts': data}, {}, CONTRACT, mode='test')
    (out / 'partial.tsv').write_text('incomplete\n')
    finished = provenance.finish_run(out, record, error=RuntimeError('deliberate probe'))
    assert finished['status'] == finished['output']['status'] == 'failed'
    assert finished['error'] == dict(type='RuntimeError', message='deliberate probe')
    assert finished['output']['files']['partial.tsv']['sha256'] == hashlib.sha256(b'incomplete\n').hexdigest()


@pytest.mark.parametrize('children', [['child-output-content-id'], {'cpu': 'child-output-content-id'}])
def test_parent_output_binds_child_evidence(tmp_path, children):
    out = tmp_path / 'parent'
    out.mkdir()
    record = provenance.start_run(out, {}, {}, CONTRACT, mode='test')
    record['children'] = children
    finished = provenance.finish_run(out, record, backend='cpu')
    assert finished['output']['children'] == record['children']
    assert finished['output']['id'] == expected_id(finished['output'])


def minimal_inputs(tmp_path):
    counts, cells, positions = (tmp_path / x for x in ['counts.h5', 'cells.tsv', 'positions.tsv'])
    with h5py.File(counts, 'w') as h:
        h.create_dataset('genes', data=np.asarray(['G'], object), dtype=h5py.string_dtype())
        h.create_dataset('cell_ids', data=np.asarray(['c'], object), dtype=h5py.string_dtype())
        h['X/data'], h['X/indices'], h['X/indptr'] = [1], [0], [0, 1]
    pd.DataFrame(dict(cell_id=['c'], specimen=['S'], patient=['P'], role=['query'], ref_group=['NA'],
                      ref_split=['NA'])).to_csv(cells, sep='\t', index=False)
    pd.DataFrame(dict(gene=['G'], chromosome=['chr1'], start=[1], end=[2])).to_csv(positions, sep='\t', index=False)
    return counts, cells, positions


def test_pipeline_existing_output_is_unchanged(tmp_path):
    from haplodose.pipeline import run
    inputs = minimal_inputs(tmp_path)
    out = tmp_path / 'out'
    out.mkdir()
    sentinel = out / 'previous_result.txt'
    sentinel.write_bytes(b'preserve existing result\n')
    with pytest.raises(FileExistsError):
        run(*inputs, out, **CONTRACT)
    assert sentinel.read_bytes() == b'preserve existing result\n'
    assert list(out.iterdir()) == [sentinel]


def test_pipeline_failure_preserves_exception_and_failure_manifest(tmp_path, monkeypatch):
    import haplodose.pipeline as pipeline
    inputs = minimal_inputs(tmp_path)
    out = tmp_path / 'out'

    def deliberate_failure(*args, **kwargs):
        (out / 'partial.tsv').write_text('partial\n')
        raise RuntimeError('execution failed at probe')

    monkeypatch.setattr(pipeline, '_run_single', deliberate_failure)
    with pytest.raises(RuntimeError, match='execution failed at probe'):
        pipeline.run(*inputs, out, **CONTRACT)
    manifest = json.loads((out / 'run_manifest.json').read_text())
    assert manifest['status'] == manifest['output']['status'] == 'failed'
    assert manifest['error']['message'] == 'execution failed at probe'
    assert 'partial.tsv' in manifest['output']['files']
