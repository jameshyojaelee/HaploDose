"""Do GEO line labels match DepMap? Spearman of each labelled line's pseudobulk (log CPM, all genes) with DepMap
bulk expression of HCT116, SW480 and RKO, after centring each gene across the three lines (identity, not tissue).
Usage: label_check.py INPUT_DIR"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from bulk_ceiling import expression  # noqa: E402

inp = Path(sys.argv[1])
with h5py.File(inp / 'counts.h5') as h:
    genes = h['genes'].asstr()[:]
    X = sparse.csr_matrix((h['X/data'][:], h['X/indices'][:], h['X/indptr'][:]), shape=(len(h['cell_ids']), len(genes)))
cells = pd.read_csv(inp / 'cells.tsv', sep='\t')
pb = pd.DataFrame({l: np.asarray(X[cells.truth_line.eq(l).to_numpy()].sum(0)).ravel() for l in ['HCT116', 'SW480', 'RKO']},
                  index=genes)
pb = np.log1p(pb / pb.sum() * 1e6)
ex = expression()
common = pb.index.intersection(ex.index)
pb, ex = pb.loc[common], ex.loc[common, ['HCT116', 'SW480', 'RKO']]
keep = (pb.mean(axis=1) > 1) & (ex.mean(axis=1) > 1)
pc, ec = pb[keep].sub(pb[keep].mean(axis=1), axis=0), ex[keep].sub(ex[keep].mean(axis=1), axis=0)
print(f'genes {keep.sum()}; rows = GEO label, columns = DepMap line')
print(pd.DataFrame({d: {g: spearmanr(pc[g], ec[d]).statistic for g in pc} for d in ec}).round(2).to_string())
