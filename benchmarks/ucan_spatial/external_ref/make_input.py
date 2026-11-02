"""Combined Haplodose input: one U-CAN section's existing input plus the external normal-epithelium tiles.

Usage: make_input.py SECTION_INPUT_DIR REFERENCE_DIR OUTDIR
Writes OUTDIR/{counts.h5, cells.tsv, positions.tsv (copy), input.json}. The section's rows come first,
unchanged; reference genes must equal the section's gene list (same panel, same order).
"""
import json
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


def read_csr(path):
    with h5py.File(path, 'r') as h:
        genes = pd.Index(h['genes'].asstr()[:])
        cells = h['cell_ids'].asstr()[:]
        X = sparse.csr_matrix((h['X/data'][:], h['X/indices'][:], h['X/indptr'][:]), shape=(len(cells), len(genes)))
    return genes, cells, X


src, ref, out = (Path(p) for p in sys.argv[1:4])
out.mkdir(parents=True, exist_ok=False)
g_s, c_s, X_s = read_csr(src / 'counts.h5')
g_r, c_r, X_r = read_csr(ref / 'counts.h5')
# Gene names repeat in the panel, so alignment by name is ambiguous: the reference was written in the
# Space Ranger feature order (aligned by Ensembl id) and must match the section's gene list exactly.
if not np.array_equal(g_s.to_numpy(), g_r.to_numpy()):
    raise ValueError('section and reference gene lists differ')
X = sparse.vstack([X_s, X_r.astype(X_s.dtype)]).tocsr()
X.sort_indices()
m_s = pd.read_csv(src / 'cells.tsv', sep='\t', dtype=str, keep_default_na=False)
m_r = pd.read_csv(ref / 'cells.tsv', sep='\t', dtype=str, keep_default_na=False)
if not np.array_equal(m_s.cell_id, c_s) or not np.array_equal(m_r.cell_id, c_r):
    raise ValueError('cells.tsv order differs from counts.h5')
if list(m_r.columns) != list(m_s.columns):
    raise ValueError(f'column mismatch {list(m_s.columns)} vs {list(m_r.columns)}')
meta = pd.concat([m_s, m_r], ignore_index=True)
if meta.cell_id.duplicated().any():
    raise ValueError('duplicated cell ids')
if not np.array_equal(np.asarray(X.sum(axis=1)).ravel()[len(c_s):], m_r.total_umi_included.astype(int)):
    raise ValueError('reference UMI totals changed in alignment')
with h5py.File(out / 'counts.h5', 'w') as h:
    h.create_dataset('genes', data=g_s.to_numpy().astype(h5py.string_dtype()))
    h.create_dataset('cell_ids', data=meta.cell_id.to_numpy(dtype=h5py.string_dtype()))
    grp = h.create_group('X')
    for key, values in [('data', X.data), ('indices', X.indices), ('indptr', X.indptr)]:
        grp.create_dataset(key, data=values, compression='gzip')
meta.to_csv(out / 'cells.tsv', sep='\t', index=False)
shutil.copyfile(src / 'positions.tsv', out / 'positions.tsv')
(out / 'input.json').write_text(json.dumps(dict(
    section_input=str(src), reference=str(ref), n_section_rows=len(c_s), n_reference_rows=len(c_r),
    reference_groups=m_r.groupby(['ref_group', 'ref_split']).size().rename('n').reset_index().to_dict('records')),
    indent=2) + '\n')
print(meta.groupby(['role', 'ref_group', 'ref_split']).size().to_string())
