#!/usr/bin/env python3
"""Claim 4 SPECTRUM: HaploDose input contract per eligible patient, from the CellRanger runs and the authors' labels.

Cells and roles (proposal in dev/spectrum_run_plan_20261103_0750/plan.md; written into docs/STATUS.md before any
tool runs). Every cell comes from the authors' GEO cell table (GSE180661_GEO_cells.tsv.gz), so all tools see the
same cells and no tool output decides who is reference or query:
- query: cell_type == Ovarian.cancer.cell, up to --max-query per sample (seeded draw); specimen = sample, patient.
- reference: the patient's other labelled cell types except 'Other' (T, B, plasma, myeloid, dendritic, mast,
  endothelial cells, fibroblasts), pooled over the patient's samples, up to --ref-per-type per cell type (seeded
  draw); ref_group = cell type; ref_split fit/holdout 50/50 (seeded). A reference cell keeps its own sample as
  specimen, so held-out reference cells sit inside each specimen (tools report their malignant rate).
- cell_type 'Other' and unsampled cells are left out.
Counts are read from CellRanger raw_feature_bc_matrix.h5 (GEO barcode + '-1'), not from CellRanger's own cell
calls, so the cell set is the authors'. GEO barcodes missing from the raw matrix are counted in qc.tsv and dropped.
Genes: CellRanger feature names (Gene Expression only; duplicated symbols keep the first) that have a position in
the GRCh38-2020-A positions.tsv used for Andor (same reference build as the CellRanger runs).
No expression value is summarised or inspected here beyond the QC counts in qc.tsv (cells matched, cells per role).
Outputs: OUT/<patient>/{counts.h5, cells.tsv, positions.tsv, barcodes/<sample>.tsv}, OUT/qc.tsv, OUT/params.json.
barcodes/<sample>.tsv is the BAM barcode list (pileup input) of that sample's cells in cells.tsv.
Usage: prepare_inputs.py SAMPLES_TSV GEO_CELLS CR_ROOT OUT [--max-query N] [--ref-per-type N] [--patients P ...]
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

POSITIONS = '/gpfs/commons/home/jameslee/CRC/results/cnv_caller_benchmark/YCC230/input/positions.tsv'
SEED = 20261103
CANCER, EXCLUDED_TYPES = 'Ovarian.cancer.cell', {'Other'}


def read_raw(h5_path, barcodes):
    """Cells x genes CSR for the given barcodes (in that order) from a CellRanger raw matrix; missing -> None."""
    with h5py.File(h5_path, 'r') as h:
        m = h['matrix']
        bcs = pd.Index(m['barcodes'][:].astype(str))
        ft = m['features/feature_type'][:].astype(str)
        names = m['features/name'][:].astype(str)
        X = sparse.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=tuple(m['shape'][:]))
    if not bcs.is_unique or not pd.Index(barcodes).is_unique:
        raise ValueError('source and requested cell barcodes must be unique')
    if X.shape != (len(names), len(bcs)) or len(ft) != len(names):
        raise ValueError('raw matrix dimensions disagree with feature/cell identifiers')
    j = bcs.get_indexer(barcodes)
    found = j >= 0
    gex = np.flatnonzero(ft == 'Gene Expression')
    return X[gex][:, j[found]].T.tocsr().astype(np.float32), names[gex], found


def align_gene_indices(names, genes):
    """First source occurrence per symbol; reject missing genes before indexing."""
    names, genes = pd.Index(names), pd.Index(genes)
    if not genes.is_unique or genes.hasnans or names.hasnans:
        raise ValueError('target genes must be unique and gene identifiers must be nonmissing')
    first = ~names.duplicated()
    match = names[first].get_indexer(genes)
    if (match < 0).any():
        raise ValueError(f'gene sets differ between samples: {list(genes[match < 0][:5])}')
    return np.flatnonzero(first)[match]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('samples')
    ap.add_argument('geo_cells')
    ap.add_argument('cr_root')
    ap.add_argument('out')
    ap.add_argument('--max-query', type=int, default=3000)
    ap.add_argument('--ref-per-type', type=int, default=300)
    ap.add_argument('--patients', nargs='*')
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=False)
    samples = pd.read_csv(a.samples, sep='\t')
    if samples.sample_title.duplicated().any():
        raise ValueError('sample manifest has duplicated sample identifiers')
    if a.patients:
        samples = samples[samples.patient.isin(a.patients)]
    geo = pd.read_csv(a.geo_cells, sep='\t', usecols=['cell_id', 'sample', 'cell_type', 'patient_id'])
    if geo.cell_id.isna().any() or geo.cell_id.duplicated().any():
        raise ValueError('GEO cell identifiers must be nonmissing and unique')
    geo = geo[geo['sample'].isin(samples.sample_title) & ~geo.cell_type.isin(EXCLUDED_TYPES)]
    expected_patient = geo['sample'].map(samples.set_index('sample_title').patient)
    if not geo.patient_id.eq(expected_patient).all():
        raise ValueError('GEO patient identifiers disagree with the RNA sample manifest')
    pos = pd.read_csv(POSITIONS, sep='\t')
    if pos.gene.isna().any() or pos.gene.duplicated().any():
        raise ValueError('position gene identifiers must be nonmissing and unique')
    qc = []
    for patient, ps in samples.groupby('patient', sort=True):
        rng = np.random.default_rng([SEED, int(patient.rsplit('-', 1)[1])])   # per patient: a subset run draws the same cells
        g = geo[geo.patient_id.eq(patient)].copy()
        g['barcode'] = [c[len(s) + 1:] + '-1' for c, s in zip(g.cell_id, g['sample'])]
        assert all(c == f'{s}_{b[:-2]}' for c, s, b in zip(g.cell_id, g['sample'], g.barcode))
        # Draw cells (metadata only): query per sample, reference per cell type pooled over the patient's samples.
        q = g[g.cell_type.eq(CANCER)]
        q = pd.concat(x.iloc[np.sort(rng.choice(len(x), min(len(x), a.max_query), replace=False))]
                      for _, x in q.groupby('sample', sort=True))
        r = g[g.cell_type.ne(CANCER)]
        r = pd.concat(x.iloc[np.sort(rng.choice(len(x), min(len(x), a.ref_per_type), replace=False))]
                      for _, x in r.groupby('cell_type', sort=True))
        q = q.assign(role='query', ref_group='', ref_split='')
        r = r.assign(role='reference', ref_group=r.cell_type, ref_split=np.where(rng.random(len(r)) < 0.5, 'fit', 'holdout'))
        cells = pd.concat([q, r]).rename(columns={'sample': 'specimen'})
        cells['patient'] = patient
        parts, metas, genes = [], [], None
        for s, c in cells.groupby('specimen', sort=True):
            h5 = Path(a.cr_root) / s / 'outs/raw_feature_bc_matrix.h5'
            X, names, found = read_raw(h5, c.barcode.to_numpy())
            first = ~pd.Index(names).duplicated()
            if genes is None:
                genes = pos.gene[pos.gene.isin(names[first])].to_numpy()
            idx = align_gene_indices(names, genes)
            parts.append(X[:, idx])
            metas.append(c[found])
            qc.append(dict(patient=patient, sample=s, geo_cells=len(c), matched=int(found.sum()),
                           query=int(c.role.eq('query')[found].sum()), reference=int(c.role.eq('reference')[found].sum()),
                           source_features=len(names), duplicate_source_symbols=int((~first).sum()),
                           gene_symbol_policy='first source occurrence; target position genes must be unique',
                           raw_h5_sha256=sha256(h5)))
        X = sparse.vstack(parts).tocsr()
        cells = pd.concat(metas, ignore_index=True)
        d = out / patient
        (d / 'barcodes').mkdir(parents=True)
        with h5py.File(d / 'counts.h5', 'w') as h:
            h.create_dataset('X/data', data=X.data.astype(np.int32), compression='gzip')
            h.create_dataset('X/indices', data=X.indices.astype(np.int32), compression='gzip')
            h.create_dataset('X/indptr', data=X.indptr.astype(np.int64))
            h.create_dataset('genes', data=genes.astype(object), dtype=h5py.string_dtype())
            h.create_dataset('cell_ids', data=cells.cell_id.to_numpy(dtype=object), dtype=h5py.string_dtype())
        cells[['cell_id', 'specimen', 'patient', 'role', 'ref_group', 'ref_split', 'cell_type', 'barcode']].to_csv(
            d / 'cells.tsv', sep='\t', index=False)
        pos[pos.gene.isin(genes)].to_csv(d / 'positions.tsv', sep='\t', index=False)
        for s, c in cells.groupby('specimen'):
            c.barcode.to_csv(d / 'barcodes' / f'{s}.tsv', index=False, header=False)
        print(patient, len(cells), 'cells', len(genes), 'genes', flush=True)
    pd.DataFrame(qc).to_csv(out / 'qc.tsv', sep='\t', index=False)
    (out / 'params.json').write_text(json.dumps(dict(vars(a), seed=SEED, positions=POSITIONS, cancer_label=CANCER,
                                                     excluded_types=sorted(EXCLUDED_TYPES),
                                                     geo_cells_sha256=sha256(a.geo_cells)), indent=2) + '\n')


if __name__ == '__main__':
    main()
