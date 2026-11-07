"""Reshape Kinker 2020 development data (results/kinker_input_v2_20261101_1551) into the Schmid et al. 2025 input
layout, for an end-to-end plumbing test of run_haplodose.py and lib_haplodose.R. Not Schmid data.

One pool, rotation r0: the query line's cells (group = DATASET) and that pool's reference-role lines (group = line
name; all listed in ref_groups.txt). Writes into OUT:
  input_<DATASET>/counts.mtx, genes.txt, cells.txt  (count_matrix.txt is then written by write_count_matrix.R,
                                                    with R write.table exactly as their preprocessing scripts do)
  input_<DATASET>/sample_annotation.txt, ref_groups.txt
  <DATASET>_WGS_groundtruth/wgs_results_formated.csv
Truth, in their read_wgs_results format (V1 index, chr without prefix, start, end, loss_wgs, base_wgs, gain_wgs):
DepMap OmicsCNGene of the query line minus the UMI-weighted mean over the reference lines (as in
benchmarks/kinker/evaluate.py), per gene; genes placed by hg38_gencode_v27.txt; 1 Mb tiles; per tile the share of
genes with difference < -0.15 (loss), > 0.15 (gain), or between (base); tiles with >= 3 genes. These thresholds are
test settings, not calibrated: the truth is relative to the reference mix, not absolute copy number.
Usage: build_kinker_input.py OUT POOL QUERY_LINE DATASET
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import io, sparse

KIN = Path('/gpfs/commons/home/jameslee/Haplodose/results/kinker_input_v2_20261101_1551')
GENCODE = Path('/gpfs/commons/home/jameslee/Haplodose/external/annotations/hg38_gencode_v27.txt')
TILE, THR, MIN_GENES = 1_000_000, 0.15, 3


def main(out, pool, query_line, dataset):
    out = Path(out)
    inp = out / f'input_{dataset}'
    inp.mkdir(parents=True)
    cells = pd.read_csv(KIN / 'cells_r0.tsv', sep='\t', dtype={'cell_id': str})
    c = cells[cells.dataset.eq(pool) & (cells.role.eq('reference') | (cells.role.eq('query') & cells.line.eq(query_line)))]
    c = c.assign(group=np.where(c.role.eq('query'), dataset, c.line))
    with h5py.File(KIN / 'counts.h5') as h:
        ids = pd.Index(h['cell_ids'][:].astype(str))
        genes = h['genes'][:].astype(str)
        X = sparse.csr_matrix((h['X/data'][:], h['X/indices'][:], h['X/indptr'][:]), shape=(len(ids), len(genes)))
    rows = ids.get_indexer(c.cell_id)
    assert (rows >= 0).all()
    Xs = X[rows].T.tocoo()                        # genes x cells
    io.mmwrite(inp / 'counts.mtx', Xs, field='integer')
    pd.Series(genes).to_csv(inp / 'genes.txt', index=False, header=False)
    c.cell_id.to_csv(inp / 'cells.txt', index=False, header=False)
    c[['cell_id', 'group']].to_csv(inp / 'sample_annotation.txt', sep='\t', index=False, header=False)
    refs = sorted(c.loc[c.role.eq('reference'), 'line'].unique())
    pd.DataFrame({'ref_groups': refs}).to_csv(inp / 'ref_groups.txt', sep='\t', index=False)

    cn = pd.read_csv(KIN / 'depmap_cn.tsv.gz', sep='\t', index_col=0)
    w = c[c.role.eq('reference')].groupby('line').umi.sum()
    truth = cn[query_line] - cn[w.index].to_numpy() @ (w / w.sum()).to_numpy()
    pos = pd.read_csv(GENCODE, sep='\t', header=None, names=['gene', 'chr', 'start', 'end', 'x'])
    pos = pos[pos.chr.isin([f'chr{i}' for i in list(range(1, 23)) + ['X']])].drop_duplicates('gene')
    g = pos.merge(truth.rename('d').dropna(), left_on='gene', right_index=True)
    g = g[g.gene.isin(genes)]
    g['tile'] = g.start // TILE
    t = g.groupby(['chr', 'tile']).agg(n=('d', 'size'), loss_wgs=('d', lambda d: (d < -THR).mean()),
                                       gain_wgs=('d', lambda d: (d > THR).mean())).reset_index()
    t = t[t.n >= MIN_GENES]
    t['base_wgs'] = 1 - t.loss_wgs - t.gain_wgs
    t['start'] = t.tile * TILE + 1
    t['end'] = (t.tile + 1) * TILE
    t['chr'] = t.chr.str.replace('chr', '', regex=False)
    order = {str(i): i for i in range(1, 23)} | {'X': 23}
    t = t.sort_values(['chr', 'start'], key=lambda s: s.map(order) if s.name == 'chr' else s)
    wgs = out / f'{dataset}_WGS_groundtruth'
    wgs.mkdir()
    t[['chr', 'start', 'end', 'loss_wgs', 'base_wgs', 'gain_wgs']].reset_index(drop=True).to_csv(wgs / 'wgs_results_formated.csv')
    print(f'cells {len(c)} (query {int(c.role.eq("query").sum())}, reference {int(c.role.eq("reference").sum())} in '
          f'{len(refs)} lines); genes {len(genes)}; truth tiles {len(t)}; '
          f'tiles with loss share > 0.5: {int((t.loss_wgs > 0.5).sum())}, gain share > 0.5: {int((t.gain_wgs > 0.5).sum())}')


if __name__ == '__main__':
    main(*sys.argv[1:5])
