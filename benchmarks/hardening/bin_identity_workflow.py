"""Public-workflow replay of equivalent named observations across old and corrected wheels."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np
import pandas as pd

from haplodose.io import CountSource
from haplodose.provenance import start_run, finish_run, file_record
from haplodose.simulate import write_dataset


RUN = r'''
import hashlib,json,sys
from pathlib import Path
import haplodose
from haplodose.pipeline import run
import haplodose._core as native
counts,cells,positions,out,artifact,source=sys.argv[1:]
artifact=json.loads(Path(artifact).read_text())
assert hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest()==artifact['native']['sha256']
assert hashlib.sha256(Path(artifact['wheel']['path']).read_bytes()).hexdigest()==artifact['wheel']['sha256']
prefix=Path(sys.executable).parent.parent
assert Path(haplodose.__file__).resolve().is_relative_to(prefix)
manifest=json.loads((Path(source)/'source_manifest.json').read_text())
for name,digest in manifest.items():
    if name.startswith('python/haplodose/') and name.endswith(('.py','.tsv')):
        actual=Path(haplodose.__file__).parent/name.removeprefix('python/haplodose/')
        assert hashlib.sha256(actual.read_bytes()).hexdigest()==digest
qc=run(counts,cells,positions,out,chunk=97,min_gene_mean=0.,device='cpu',
    counts_kind='raw',genome_build='hg38',coordinate_system='one-based-closed',gene_namespace='symbol')
assert qc['device']==qc['backend']=='cpu'
print(json.dumps(dict(wheel=artifact['wheel'],package=str(haplodose.__file__),backend=qc['backend'])))
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    old_task = args.task.parent/'spatial_identity_20261109'
    old_artifact = old_task/'validation_23030652/installed_artifact.json'
    current_artifact = Path(os.environ['HAPLODOSE_VALIDATED_ARTIFACT'])
    inputs = dict(script=Path(__file__), design=args.task/'WORKFLOW_DESIGN.md',
        old_artifact=old_artifact, candidate_artifact=current_artifact,
        old_source=args.task/'probe_source_v2/source_manifest.json',
        candidate_source=args.task/'candidate_v1/source_manifest.json')
    run = start_run(args.out, inputs, dict(seed=20261031, min_gene_mean=0., column_order='original_and_reverse'),
        dict(scope='Four simulated identity replays, no biological or performance calibration'), mode='bin_identity_workflow')
    try:
        data = write_dataset(args.out/'data', seed=20261031)
        positions = pd.read_csv(data/'positions.tsv',sep='\t')
        assert (positions.gene=='G00009').sum()==(positions.gene=='G00010').sum()==1
        positions.loc[positions.gene=='G00010','start'] = positions.loc[positions.gene=='G00009','start'].iloc[0]
        positions.to_csv(args.out/'tied_positions.tsv',sep='\t',index=False)
        with CountSource(data/'counts.h5',gene_namespace='symbol') as source:
            genes = source.genes.copy(); ids = source.cell_ids.copy(); counts = source.read(0,source.n_cells)
        order = np.arange(len(genes))[::-1]; reordered = counts[:,order].tocsr(); reordered.sort_indices()
        assert np.array_equal(reordered[:,order].toarray(),counts.toarray())
        with h5py.File(args.out/'reversed_counts.h5','x') as handle:
            handle.create_dataset('X/data',data=reordered.data.astype(np.int32),compression='gzip')
            handle.create_dataset('X/indices',data=reordered.indices.astype(np.int32),compression='gzip')
            handle.create_dataset('X/indptr',data=reordered.indptr.astype(np.int64))
            handle.create_dataset('genes',data=genes[order].astype(object),dtype=h5py.string_dtype())
            handle.create_dataset('cell_ids',data=ids.astype(object),dtype=h5py.string_dtype())
        (args.out/'input_identities.json').write_text(json.dumps({str(path):file_record(path) for path in
            (data/'counts.h5',data/'cells.tsv',args.out/'tied_positions.tsv',args.out/'reversed_counts.h5')},indent=2)+'\n')
        for arm, python, artifact, source in (
            ('old',old_task/'environment/bin/python',old_artifact,args.task/'probe_source_v2'),
            ('corrected',Path(sys.executable),current_artifact,args.task/'candidate_v1')):
            env = {**os.environ,'LD_LIBRARY_PATH':str(python.parent.parent/'lib'),'JAX_PLATFORMS':'cpu'}
            for layout, matrix in (('original',data/'counts.h5'),('reversed',args.out/'reversed_counts.h5')):
                with (args.out/f'{arm}_{layout}.stdout').open('w') as stdout,(args.out/f'{arm}_{layout}.stderr').open('w') as stderr:
                    subprocess.run([str(python),'-I','-c',RUN,str(matrix),str(data/'cells.tsv'),
                        str(args.out/'tied_positions.tsv'),str(args.out/f'{arm}_{layout}'),str(artifact),str(source)],
                        env=env,stdout=stdout,stderr=stderr,check=True,timeout=1800)
        before_a, before_b = args.out/'old_original',args.out/'old_reversed'
        after_a, after_b = args.out/'corrected_original',args.out/'corrected_reversed'
        assert not np.array_equal(np.load(before_a/'bin_counts.npy'),np.load(before_b/'bin_counts.npy'))
        arrays = ['bin_counts.npy','posteriors.npy']
        for name in arrays:
            first, second = np.load(after_a/name), np.load(after_b/name)
            assert first.shape == second.shape and first.dtype == second.dtype, name
            missing = np.isnan(first)
            assert np.array_equal(missing,np.isnan(second)), name
            assert np.array_equal(first[~missing],second[~missing]), name
        tables = ['bins.tsv','cells.tsv.gz','profiles.tsv.gz','segments.tsv','differences.tsv','clones.tsv','arms.tsv','reference_check.tsv']
        for name in tables:
            pd.testing.assert_frame_equal(pd.read_csv(after_a/name,sep='\t'),pd.read_csv(after_b/name,sep='\t'),check_exact=True)
        first = pd.read_csv(after_a/'genes.tsv.gz',sep='\t').drop(columns='column').set_index('gene').sort_index()
        second = pd.read_csv(after_b/'genes.tsv.gz',sep='\t').drop(columns='column').set_index('gene').sort_index()
        pd.testing.assert_frame_equal(first,second,check_exact=True)
        (args.out/'report.json').write_text(json.dumps(dict(fits=4, biological_specimens=0,
            original_wheel_bin_counts_differ=True, corrected_equivalent_arrays_exact=arrays,
            corrected_equivalent_tables_exact=tables, gene_metadata_identity_aligned_exact=True,
            actual_backends_inside_inference='cpu', source_and_native_identities_verified=True,
            scientific_thresholds_changed=False, recovery_or_error_control_claim=False),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error)
        raise


if __name__=='__main__': main()
