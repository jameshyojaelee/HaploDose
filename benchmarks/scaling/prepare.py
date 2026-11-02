"""Inputs for the compute scaling benchmark: real Kinker 2020 cells at several sizes, same files for every tool.

Source: one Kinker rotation (prepare.py output: counts.h5, cells_r0.tsv, positions.tsv). Only query and reference
cells are used (every tool analyses every cell it is given). Sizes up to the available cells are drawn without
replacement, stratified by role and reference split; larger sizes are drawn with replacement and duplicates get
new ids (<cell_id>_d<k>), so expression is real but cells repeat. Each size gets its own counts_<N>.h5 (written in
50k-cell chunks) and cells_<N>.tsv; positions.tsv is copied. Seed 20261102.
Usage: prepare.py INPUT_DIR OUT_DIR SIZE [SIZE ...]
"""
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

SEED = 20261102
CHUNK = 50_000


def read_h5(path):
    with h5py.File(path) as h:
        X = sparse.csr_matrix((h['X/data'][:], h['X/indices'][:], h['X/indptr'][:]))
        return X, h['genes'][:].astype(str), h['cell_ids'][:].astype(str)


def write_h5(path, X, rows, genes, cell_ids):
    """CSR rows X[rows] in chunks, so a 1M-cell file never holds all rows in memory at once."""
    with h5py.File(path, 'w') as h:
        h.create_dataset('genes', data=genes.astype(h5py.string_dtype()))
        h.create_dataset('cell_ids', data=np.asarray(cell_ids).astype(h5py.string_dtype()))
        g = h.create_group('X')
        data = g.create_dataset('data', shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(1 << 20,),
                                compression='gzip')
        idx = g.create_dataset('indices', shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(1 << 20,),
                               compression='gzip')
        indptr = [np.zeros(1, np.int64)]
        for s in range(0, len(rows), CHUNK):
            B = X[rows[s:s + CHUNK]]
            k = data.shape[0]
            data.resize((k + B.nnz,))
            idx.resize((k + B.nnz,))
            data[k:] = B.data.astype(np.int32)
            idx[k:] = B.indices.astype(np.int32)
            indptr.append(k + B.indptr[1:].astype(np.int64))
        g.create_dataset('indptr', data=np.concatenate(indptr))


def main(inp, out, *sizes):
    inp, out = Path(inp), Path(out)
    out.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(SEED)
    X, genes, ids = read_h5(inp / 'counts.h5')
    cells = pd.read_csv(inp / 'cells_r0.tsv', sep='\t', dtype=str, keep_default_na=False)
    cells = cells[cells.role.isin(['query', 'reference'])].reset_index(drop=True)
    row = pd.Series(np.arange(len(ids)), index=ids).loc[cells.cell_id].to_numpy()
    strata = (cells.role + '|' + cells.ref_split).to_numpy()
    shutil.copy(inp / 'positions.tsv', out / 'positions.tsv')
    log = []
    for n in map(int, sizes):
        take = []
        for s in np.unique(strata):
            members = np.flatnonzero(strata == s)
            k = int(round(n * len(members) / len(cells)))
            take.append(rng.choice(members, k, replace=n > len(cells)))
        take = np.sort(np.concatenate(take))
        sub = cells.iloc[take].copy()
        dup = sub.groupby('cell_id').cumcount().to_numpy()
        sub['cell_id'] = np.where(dup > 0, sub.cell_id + '_d' + dup.astype(str), sub.cell_id)
        sub.to_csv(out / f'cells_{n}.tsv', sep='\t', index=False)
        write_h5(out / f'counts_{n}.h5', X, row[take], genes, sub.cell_id.to_numpy())
        log.append(dict(size=n, cells=len(sub), unique_source_cells=int(len(np.unique(take))),
                        specimens=int(sub.specimen.replace('', np.nan).nunique())))
        print(log[-1], flush=True)
    pd.DataFrame(log).to_csv(out / 'sizes.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main(*sys.argv[1:])
