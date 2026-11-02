"""Kinker 2020 (GSE157220) pooled cell-line scRNA input, with cell-line labels rebuilt and DepMap truth.

GEO holds linear CPM only and no per-cell line labels (those are on the Broad Single Cell Portal, login
needed). Steps:
1. Counts: per cell, library = 1e6 / smallest nonzero CPM; counts = CPM x library / 1e6. Every recovered
   value must be within 1e-3 of an integer, otherwise the script stops.
2. Pools: cell name suffix -<pool> (8 pools of 21-25 lines); names c<N> are the 8-line 'custom' pool.
3. Line per cell, within its pool's listed lines (GSE157220_Pool_composition.xlsx, exported to TSV): Pearson
   correlation of the cell's centred log1p(CP10k) with each line's centred DepMap 24Q4 expression over the
   2,000 genes most variable between the pool's lines; argmax, margin = best - second. Checked against the
   published post-QC cells per line.
4. Design, two rotations per pool. Eligible lines (see eligible()) are split in half at random (seed); in rotation r,
   half r is the reference (one ref_group per pool, 70% fit / 30% holdout) and every line of the other half
   is its own single-line specimen. Each line is a query exactly once. dataset = pool, so a pool's specimens
   use only that pool's reference lines.
This is expression from cancer lines against other cancer lines, not tumour against normal: the reference
carries its own CN, which the truth (line CN minus the reference lines' CN) accounts for.
Usage: prepare.py POOL_TSV OUT
"""
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

SRC = Path('/gpfs/commons/home/jameslee/Haplodose/data/external/kinker2020_gse157220/GSE157220_CPM_data.txt.gz')
DEPMAP = Path('/gpfs/commons/home/jameslee/CRC/data/external/depmap_24q4')
POSITIONS = Path('/gpfs/commons/home/jameslee/CRC/results/cnv_caller_benchmark/SMC_full/input/positions.tsv')
SEED = 20261101
N_GENES = 2000
MIN_MARGIN = 0.05       # cells below this assignment margin are left out of the design


def read_counts():
    blocks, genes = [], []
    reader = pd.read_csv(SRC, sep='\t', index_col=0, chunksize=1000, engine='c')
    for chunk in reader:
        genes.extend(chunk.index.astype(str))
        blocks.append(sparse.csr_matrix(chunk.to_numpy()))
        cells = chunk.columns.astype(str)
    X = sparse.vstack(blocks).T.tocsr()                     # cells x genes, CPM
    mn = np.minimum.reduceat(X.data, X.indptr[:-1])        # every cell has nonzero genes
    lib = 1e6 / mn
    data = X.data * np.repeat(lib, np.diff(X.indptr)) / 1e6
    dev = float(np.abs(data - np.round(data)).max())
    if dev > 1e-3:
        raise ValueError(f'CPM does not map back to integer counts (max deviation {dev})')
    X.data = np.round(data)
    return X.astype(np.int32), np.array(genes), np.array(cells), lib, dev


def depmap_rows(path, ids):
    header = subprocess.run(['head', '-1', str(path)], capture_output=True, text=True, check=True).stdout
    rows = subprocess.run(['grep', '-E', '^(' + '|'.join(ids) + '),', str(path)],
                          capture_output=True, text=True, check=True).stdout.splitlines()
    cols = [c.split(' (')[0] for c in header.rstrip('\n').split(',')[1:]]
    d = pd.DataFrame([[float(x) if x else np.nan for x in r.split(',')[1:]] for r in rows], columns=cols,
                     index=[r.split(',')[0] for r in rows])
    return d.loc[:, ~d.columns.duplicated()]


def assign(X, genes, lib, cells_pool, pools, expr):
    out = pd.Series('', index=range(X.shape[0]), dtype=object)
    margin = np.zeros(X.shape[0])
    gidx = pd.Index(genes)
    for p, g in pools.groupby('pool'):
        lines = [l for l in g.line if l in expr.index]
        rows = np.flatnonzero(cells_pool == p)
        e = expr.loc[lines]
        e = e.loc[:, e.columns.isin(gidx)]
        top = e.var(axis=0).sort_values(ascending=False).index[:N_GENES]
        E = e[top].to_numpy()
        E = E - E.mean(axis=0)
        E /= np.linalg.norm(E, axis=1, keepdims=True)
        C = X[rows][:, gidx.get_indexer(top)].toarray() / lib[rows, None] * 1e4
        C = np.log1p(C)
        C -= C.mean(axis=0)
        C /= np.maximum(np.linalg.norm(C, axis=1, keepdims=True), 1e-9)
        R = C @ E.T                                         # cells x lines
        o = np.sort(R, axis=1)
        out.iloc[rows] = np.array(lines)[R.argmax(axis=1)]
        margin[rows] = o[:, -1] - o[:, -2]
    return out.to_numpy(), margin


def eligible(check):
    """Lines usable in the design: DepMap CN and expression, assigned cells within 10% of the published
    count, and not the custom pool (5 of its 8 lines have no DepMap expression, so their cells were
    assigned to the other 3). A line with no expression profile sends its cells to its nearest lines,
    which then show an excess and are dropped too (pool 10: NCIH2077 -> NCIH1048, RD)."""
    ratio = check.assigned_cells / check.published_cells
    return check.has_cn & check.has_expr & ratio.between(0.9, 1.1) & check.pool.ne('custom')


def write_design(cells_t, check, out):
    rng = np.random.default_rng(SEED)
    for r in (0, 1):
        c = cells_t.copy()
        c['role'], c['ref_group'], c['ref_split'], c['specimen'], c['design'] = 'other', '', '', '', 'unused'
        c['dataset'] = 'pool' + c.pool
        usable = c.margin >= MIN_MARGIN
        for p, g in check[eligible(check)].groupby('pool'):
            names = sorted(g.line)
            perm = np.random.default_rng(SEED + int(p) if p.isdigit() else SEED + 999).permutation(len(names))
            half = {names[i] for i in perm[:len(names) // 2]}
            ref_lines = half if r == 0 else set(names) - half
            ref = usable & c.pool.eq(p) & c.line.isin(ref_lines)
            q = usable & c.pool.eq(p) & c.line.isin(set(names) - ref_lines)
            c.loc[ref, 'role'], c.loc[ref, 'ref_group'], c.loc[ref, 'design'] = 'reference', 'pool_ref', 'reference'
            c.loc[ref, 'ref_split'] = np.where(rng.random(ref.sum()) < 0.7, 'fit', 'holdout')
            c.loc[ref, 'specimen'] = f'pool{p}_ref'
            c.loc[q, 'role'], c.loc[q, 'design'] = 'query', 'single_line'
            c.loc[q, 'specimen'] = c.loc[q, 'line']
        c['patient'] = c.specimen
        c['truth_line'] = c.line
        c.to_csv(out / f'cells_r{r}.tsv', sep='\t', index=False)
        print(f'rotation {r}:', c.groupby('design').size().to_dict(), flush=True)


def main(pool_tsv, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    pools = pd.read_csv(pool_tsv, sep='\t', dtype=str).rename(columns={'Pool_ID': 'pool', 'Cell line name': 'line'})
    pools['published_cells'] = pools.Cells_Analyzed_Post_QC.astype(int)
    model = pd.read_csv(DEPMAP / 'Model.csv', usecols=['ModelID', 'CCLEName'], dtype=str).dropna()
    pools = pools.merge(model, left_on='line', right_on='CCLEName', how='left')
    X, genes, cells, lib, dev = read_counts()
    print(f'cells {len(cells)}, genes {len(genes)}, max deviation from integer {dev:.2e}, median library {np.median(lib):.0f}', flush=True)
    pool = np.array([c.rsplit('-', 1)[1] if '-' in c else 'custom' for c in cells])
    ids = pools.ModelID.dropna().tolist()
    expr = depmap_rows(DEPMAP / 'OmicsExpressionProteinCodingGenesTPMLogp1.csv', ids)
    cn = depmap_rows(DEPMAP / 'OmicsCNGene.csv', ids)
    id2line = pools.dropna(subset=['ModelID']).set_index('ModelID').line
    expr.index, cn.index = expr.index.map(id2line), cn.index.map(id2line)
    line, margin = assign(X, genes, lib, pool, pools, expr)
    umi = np.asarray(X.sum(axis=1)).ravel()
    cells_t = pd.DataFrame(dict(cell_id=cells, pool=pool, line=line, margin=margin, library=lib, umi=umi))
    check = cells_t.groupby(['pool', 'line']).size().rename('assigned_cells').reset_index()
    check = pools.merge(check, on=['pool', 'line'], how='left').fillna({'assigned_cells': 0})
    check['has_expr'] = check.line.isin(expr.index)
    check['has_cn'] = check.line.isin(cn.index)
    check.to_csv(out / 'assignment_check.tsv', sep='\t', index=False)
    print(check[['pool', 'line', 'published_cells', 'assigned_cells', 'has_expr', 'has_cn']].to_string(), flush=True)
    print('Spearman published vs assigned:', check[['published_cells', 'assigned_cells']].corr('spearman').iloc[0, 1])

    write_design(cells_t, check, out)
    cells_t.to_csv(out / 'cells_assigned.tsv.gz', sep='\t', index=False)
    cn.T.to_csv(out / 'depmap_cn.tsv.gz', sep='\t')
    pos = pd.read_csv(POSITIONS, sep='\t')
    pos[pos.gene.isin(genes)].to_csv(out / 'positions.tsv', sep='\t', index=False)
    with h5py.File(out / 'counts.h5', 'w') as h:
        h.create_dataset('genes', data=genes.astype(h5py.string_dtype()))
        h.create_dataset('cell_ids', data=cells.astype(h5py.string_dtype()))
        g = h.create_group('X')
        for k, v in [('data', X.data.astype(np.int32)), ('indices', X.indices.astype(np.int32)), ('indptr', X.indptr)]:
            g.create_dataset(k, data=v, compression='gzip')
    print(f'genes with GRCh38 position {pos.gene.isin(genes).sum()} of {len(genes)}; median UMI/cell {np.median(umi):.0f}')


if __name__ == '__main__':
    main(*sys.argv[1:])
