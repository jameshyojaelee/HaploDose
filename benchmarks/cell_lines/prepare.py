"""Cell-line benchmark input from GSE149224 (HCT116, SW480, RKO; 5-FU 0/10/50/200 uM; 10 shared batches).

The three lines were captured together in each batch, so line differences carry no platform confound.
Design (one Haplodose run):
  - RKO cells are the reference (ref_group RKO, 70% fit / 30% holdout).
  - Batches 1-5: one specimen per batch with HCT116 + SW480 cells (role query). Truth: two clones whose
    pooled log-rate difference should follow DepMap CN(HCT116) - CN(SW480).
  - Batches 6-10: one specimen per line (HCT116 or SW480 alone) pooled over these batches. Each batch holds
    one dose per line, so a specimen mixes batches and 5-FU doses. Truth: one clone; batch and dose
    (a strong expression programme) must not split it.
Positions: GRCh38 gene coordinates from the CRC SMC prepared input (symbol match; GEO symbols are older).
Usage: prepare.py OUT
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

SRC = Path('/gpfs/commons/groups/innovation/CRC_data/sc_data/PublicData_sc/GSE149224')
POSITIONS = Path('/gpfs/commons/home/jameslee/CRC/results/cnv_caller_benchmark/SMC_full/input/positions.tsv')
SEED = 20261101


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    meta = pd.read_csv(SRC / 'GSE149224_meta.information.csv', index_col=0)
    meta.columns = ['batch', 'line', 'dose']
    X = pd.read_csv(SRC / 'GSE149224_RSH.all.counts.txt', sep=' ', index_col=0, dtype={c: np.int32 for c in meta.index})
    X = X[meta.index]                                                  # genes x cells, matrix order = meta order
    genes = X.index.astype(str).to_numpy()
    M = sparse.csr_matrix(X.to_numpy().T)                             # cells x genes
    del X
    rng = np.random.default_rng(SEED)
    batches = sorted(meta.batch.unique())
    mix, single = batches[:5], batches[5:]
    cells = pd.DataFrame({'cell_id': meta.index, 'batch': meta.batch.to_numpy(), 'truth_line': meta.line.to_numpy(),
                          'dose': meta.dose.to_numpy()})
    cells['role'], cells['ref_group'], cells['ref_split'], cells['specimen'] = 'other', '', '', ''
    rko = cells.truth_line.eq('RKO')
    cells.loc[rko, 'role'] = 'reference'
    cells.loc[rko, 'ref_group'] = 'RKO'
    cells.loc[rko, 'ref_split'] = np.where(rng.random(rko.sum()) < 0.7, 'fit', 'holdout')
    cells.loc[rko, 'specimen'] = 'RKO_ref'
    q = ~rko & cells.batch.isin(mix)
    cells.loc[q, 'specimen'] = 'mix_' + cells.loc[q, 'batch'].astype(str)
    s = ~rko & cells.batch.isin(single)
    cells.loc[s, 'specimen'] = cells.loc[s, 'truth_line'] + '_pooled'
    cells.loc[q | s, 'role'] = 'query'
    cells['patient'] = cells.specimen
    cells['design'] = np.select([q, s], ['mixture', 'single_line'], 'reference')
    cols = ['cell_id', 'specimen', 'patient', 'role', 'ref_group', 'ref_split', 'truth_line', 'dose', 'batch', 'design']
    cells[cols].to_csv(out / 'cells.tsv', sep='\t', index=False)
    pos = pd.read_csv(POSITIONS, sep='\t')
    pos[pos.gene.isin(genes)].to_csv(out / 'positions.tsv', sep='\t', index=False)
    with h5py.File(out / 'counts.h5', 'w') as h:
        h.create_dataset('genes', data=genes.astype(h5py.string_dtype()))
        h.create_dataset('cell_ids', data=cells.cell_id.to_numpy().astype(h5py.string_dtype()))
        g = h.create_group('X')
        for k, v in [('data', M.data.astype(np.int32)), ('indices', M.indices.astype(np.int32)), ('indptr', M.indptr)]:
            g.create_dataset(k, data=v, compression='gzip')
    umi = np.asarray(M.sum(axis=1)).ravel()
    print(cells.groupby(['design', 'truth_line']).size().to_string())
    print(f'genes {len(genes)}, with GRCh38 position {pos.gene.isin(genes).sum()}, median UMI/cell {np.median(umi):.0f}')


if __name__ == '__main__':
    main(sys.argv[1])
