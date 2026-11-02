"""MIX-seq (McFarland et al. 2020) held-out input in the Kinker benchmark format (docs/heldout_protocol.md).

Samples (protocol amendment 1): DMSO_24hr_expt1, DMSO_6hr_expt3, DMSO_expt10 (one per pool composition).
Cells: the authors' QC-passing singlets (cell_quality == normal); line = singlet_ID (SNP-based, CCLE name),
DepMap_ID from the same table. Counts: 10x UMI matrices (genes.tsv: Ensembl ID, symbol); genes are summed by
symbol. Lines: >= 50 QC-passing cells in the sample and DepMap 24Q4 copy number; a line in several pools is
used only in the pool where it has the most cells (elsewhere its cells are 'unused'). Design as
benchmarks/kinker/prepare.py write_design (two rotations, half the eligible lines of a pool as its reference,
70% fit / 30% holdout), seed 20261103. Outputs: counts.h5, cells_r{0,1}.tsv, positions.tsv, depmap_cn.tsv.gz,
lines.tsv (per line: pool, cells, eligibility, in_kinker for the protocol's sensitivity population).
Usage: prepare_mixseq.py EXTRACTED_DIR OUT
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import io as sio
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'kinker'))
from prepare import DEPMAP, POSITIONS, depmap_rows  # noqa: E402

SEED = 20261103
SAMPLES = {'expt1': 'DMSO_24hr_expt1', 'expt3': 'DMSO_6hr_expt3', 'expt10': 'DMSO_expt10'}
MIN_CELLS = 50
KINKER_CELLS = Path(__file__).resolve().parents[2] / 'results/kinker_input_v2_20261101_1551/cells_r0.tsv'


def read_sample(d):
    X = sio.mmread(d / 'matrix.mtx').T.tocsr()                       # cells x genes
    genes = pd.read_csv(d / 'genes.tsv', sep='\t', header=None, names=['ensembl', 'symbol'])
    bc = pd.read_csv(d / 'barcodes.tsv', sep='\t', header=None)[0].astype(str).to_numpy()
    meta = pd.read_csv(d / 'classifications.csv').set_index('barcode').reindex(bc)
    return X, genes.symbol.astype(str).to_numpy(), bc, meta


def main(src, out):
    src, out = Path(src), Path(out)
    out.mkdir(parents=True, exist_ok=False)
    mats, metas, sym_all = [], [], None
    for pool, name in SAMPLES.items():
        X, sym, bc, meta = read_sample(src / name / name)
        keep = meta.cell_quality.eq('normal').to_numpy()
        X, bc, meta = X[keep], bc[keep], meta[keep]
        u, inv = np.unique(sym, return_inverse=True)                  # sum duplicate symbols
        S = sparse.csr_matrix((np.ones(len(sym)), (np.arange(len(sym)), inv)), shape=(len(sym), len(u)))
        X = (X @ S).tocsr()
        if sym_all is None:
            sym_all = u
        if not np.array_equal(u, sym_all):
            raise ValueError(f'{name}: gene list differs from the first sample')
        metas.append(pd.DataFrame(dict(cell_id=[f'{pool}_{b}' for b in bc], pool=pool, line=meta.singlet_ID.to_numpy(),
                                       depmap_id=meta.DepMap_ID.to_numpy(), margin=1.0)))
        mats.append(X)
        print(name, X.shape, flush=True)
    X = sparse.vstack(mats).tocsr().astype(np.int32)
    cells = pd.concat(metas, ignore_index=True)
    cells['umi'] = np.asarray(X.sum(axis=1)).ravel()
    cells['library'] = cells.umi.astype(float)

    model = pd.read_csv(DEPMAP / 'Model.csv', usecols=['ModelID', 'CCLEName'], dtype=str).dropna()
    ids = sorted(set(cells.depmap_id.dropna()) & set(model.ModelID))
    cn = depmap_rows(DEPMAP / 'OmicsCNGene.csv', ids)
    id2line = model.set_index('ModelID').CCLEName
    cn.index = cn.index.map(id2line)
    n = cells.groupby(['line', 'pool']).size().rename('cells').reset_index()
    n['has_cn'] = n.line.isin(cn.index)
    n['enough'] = n.cells >= MIN_CELLS
    best = n[n.has_cn & n.enough].sort_values('cells', ascending=False).drop_duplicates('line')
    n['chosen'] = n.set_index(['line', 'pool']).index.isin(best.set_index(['line', 'pool']).index)
    kinker = set(pd.read_csv(KINKER_CELLS, sep='\t', usecols=['line']).line)
    n['in_kinker'] = n.line.isin(kinker)
    n.to_csv(out / 'lines.tsv', sep='\t', index=False)
    print(f'lines chosen {int(n.chosen.sum())}; in Kinker {int((n.chosen & n.in_kinker).sum())}', flush=True)

    chosen = set(zip(n.loc[n.chosen, 'line'], n.loc[n.chosen, 'pool']))
    rng = np.random.default_rng(SEED)
    for r in (0, 1):
        c = cells.copy()
        c['role'], c['ref_group'], c['ref_split'], c['specimen'], c['design'] = 'other', '', '', '', 'unused'
        c['dataset'] = 'pool' + c.pool
        for k, (p, g) in enumerate(n[n.chosen].groupby('pool')):
            names = sorted(g.line)
            perm = np.random.default_rng(SEED + k).permutation(len(names))
            half = {names[i] for i in perm[:len(names) // 2]}
            ref_lines = half if r == 0 else set(names) - half
            ref = c.pool.eq(p) & c.line.isin(ref_lines)
            q = c.pool.eq(p) & c.line.isin(set(names) - ref_lines)
            c.loc[ref, 'role'], c.loc[ref, 'ref_group'], c.loc[ref, 'design'] = 'reference', 'pool_ref', 'reference'
            c.loc[ref, 'ref_split'] = np.where(rng.random(ref.sum()) < 0.7, 'fit', 'holdout')
            c.loc[ref, 'specimen'] = f'pool{p}_ref'
            c.loc[q, 'role'], c.loc[q, 'design'] = 'query', 'single_line'
            c.loc[q, 'specimen'] = c.loc[q, 'line']
        c['patient'] = c.specimen
        c['truth_line'] = c.line
        assert not c[c.design.eq('single_line')].groupby('line').pool.nunique().gt(1).any()
        c.drop(columns='depmap_id').to_csv(out / f'cells_r{r}.tsv', sep='\t', index=False)
        print(f'rotation {r}:', c.groupby('design').size().to_dict(), flush=True)
    cn.T.to_csv(out / 'depmap_cn.tsv.gz', sep='\t')
    pos = pd.read_csv(POSITIONS, sep='\t')
    pos[pos.gene.isin(sym_all)].to_csv(out / 'positions.tsv', sep='\t', index=False)
    with h5py.File(out / 'counts.h5', 'w') as h:
        h.create_dataset('genes', data=sym_all.astype(h5py.string_dtype()))
        h.create_dataset('cell_ids', data=cells.cell_id.to_numpy().astype(h5py.string_dtype()))
        g = h.create_group('X')
        for k, v in [('data', X.data.astype(np.int32)), ('indices', X.indices.astype(np.int32)), ('indptr', X.indptr)]:
            g.create_dataset(k, data=v, compression='gzip')
    print(f'genes with position {pos.gene.isin(sym_all).sum()} of {len(sym_all)}; median UMI/cell {cells.umi.median():.0f}')


if __name__ == '__main__':
    main(*sys.argv[1:])
