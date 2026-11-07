"""Stream independent dataset inputs without sharing selected genomic features."""
from contextlib import ExitStack
import hashlib
from pathlib import Path

import h5py
import numpy as np

from .io import CountSource, write_table


def materialize_datasets(counts, meta, directory, *, gene_namespace, chunk):
    directory = Path(directory)
    directory.mkdir(exist_ok=False)
    output = []
    with CountSource(counts, gene_namespace=gene_namespace) as src, ExitStack() as stack:
        handles = {}
        labels = meta.dataset.to_numpy()
        for i, dataset in enumerate(sorted(meta.dataset.unique())):
            key = f'dataset_{i:03d}_{hashlib.sha256(dataset.encode()).hexdigest()[:10]}'
            chosen = meta.dataset.eq(dataset).to_numpy()
            path = directory / f'{key}.h5'
            h = stack.enter_context(h5py.File(path, 'x'))
            g = h.create_group('X')
            for name, dtype in [('data', np.float32), ('indices', np.int32)]:
                g.create_dataset(name, (0,), maxshape=(None,), dtype=dtype, chunks=True)
            g.create_dataset('indptr', data=np.array([0], np.int64), maxshape=(None,), chunks=True)
            h.create_dataset('genes', data=np.asarray(src.genes, dtype='S'))
            h.create_dataset('cell_ids', data=np.asarray(meta.loc[chosen, 'cell_id'], dtype='S'))
            cells = directory / f'{key}.tsv'
            write_table(meta.loc[chosen].reset_index(drop=True), cells)
            handles[dataset] = g
            output.append(dict(dataset=dataset, key=key, counts=path, cells=cells))
        for a, b, matrix in src.chunks(chunk):
            for dataset in np.unique(labels[a:b]):
                selected = matrix[labels[a:b] == dataset]
                g = handles[dataset]
                offset = len(g['data'])
                for name, values in [('data', selected.data), ('indices', selected.indices)]:
                    old = len(g[name])
                    g[name].resize((old + len(values),))
                    g[name][old:] = values
                old = len(g['indptr'])
                g['indptr'].resize((old + selected.shape[0],))
                g['indptr'][old:] = selected.indptr[1:] + offset
    return output
