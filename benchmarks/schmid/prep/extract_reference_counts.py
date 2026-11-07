"""Extract the fixed external reference cells (STATUS item 5, amendment f) from an h5ad as MatrixMarket.

Reads only the rows listed in the saved selection table (results/schmid_reference_cells_20261102_2119/
<DS>_reference_cells.tsv, column atlas_row is 1-based obs order), in the table's row order, which is the column
order the authors' R scripts produce (COLO320: R sample() order; MCF7: atlas order). Writes genes x cells
matrix.mtx (integer), genes.tsv (var index = gene symbol, as rownames(X) in their readH5AD call) and cells.tsv
(cell_id, group). Their R script then intersects genes and cbinds (prep_external_ref.R). Counts are checked to be
integers. Usage: extract_reference_counts.py H5AD LAYER SELECTION_TSV OUTDIR   (LAYER: X or a layers/<name> key)
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def decode(v):
    return np.array([x.decode('utf8') if isinstance(x, bytes) else str(x) for x in v])


def main(h5ad, layer, selection, outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=False)
    sel = pd.read_csv(selection, sep='\t')
    rows = sel.atlas_row.to_numpy() - 1
    with h5py.File(h5ad) as h:
        obs_idx = h['obs'].attrs.get('_index', '_index')
        obs_names = decode(h['obs'][obs_idx][:])
        assert (obs_names[rows] == sel.cell_id.astype(str).to_numpy()).all(), 'cell IDs do not match atlas rows'
        var_idx = h['var'].attrs.get('_index', '_index')
        genes = decode(h['var'][var_idx][:])
        g = h['X'] if layer == 'X' else h['layers'][layer]
        enc = g.attrs.get('encoding-type', g.attrs.get('h5sparse_format', ''))
        enc = enc.decode() if isinstance(enc, bytes) else str(enc)
        assert 'csr' in enc, f'expected CSR (cells in rows), got {enc!r}'
        indptr = g['indptr'][:]
        data, indices = g['data'], g['indices']
        parts_d, parts_i, ptr = [], [], [0]
        for r in rows:
            a, b = indptr[r], indptr[r + 1]
            parts_d.append(data[a:b])
            parts_i.append(indices[a:b])
            ptr.append(ptr[-1] + (b - a))
    d = np.concatenate(parts_d)
    assert np.all(d == np.round(d)), 'non-integer counts'
    m = sp.csr_matrix((d.astype(np.int64), np.concatenate(parts_i), np.array(ptr)), shape=(len(rows), len(genes)))
    scipy.io.mmwrite(out / 'matrix.mtx', m.T.tocoo(), field='integer')
    pd.Series(genes).to_csv(out / 'genes.tsv', sep='\t', index=False, header=False)
    sel[['cell_id', 'group']].to_csv(out / 'cells.tsv', sep='\t', index=False, header=False)
    info = dict(h5ad=str(h5ad), layer=layer, selection=str(selection), n_cells=len(rows), n_genes=len(genes),
                encoding=enc, total_counts=int(d.sum()), duplicated_gene_symbols=int(pd.Series(genes).duplicated().sum()))
    (out / 'extract.json').write_text(json.dumps(info, indent=1) + '\n')
    print(info)


if __name__ == '__main__':
    main(*sys.argv[1:5])
