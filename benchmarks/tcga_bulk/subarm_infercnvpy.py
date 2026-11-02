"""infercnvpy 0.6.1 window values per sample for subarm.py (scrna_cnv env). Same preprocessing as
benchmarks/ucan_bulk/run_infercnvpy_bulk.py (normalize_total 1e4, log1p, normal reference, window 100, step 10,
genes with mean normal CPM >= 1). Writes OUT_DIR/infercnvpy_windows.npz (Z samples x windows, samples) and
OUT_DIR/infercnvpy_gene_window.tsv (gene -> window, infercnvpy gene order).
Usage: subarm_infercnvpy.py TYPE_INPUT_DIR OUT_DIR
"""
import os
import sys
import time
from pathlib import Path

import anndata as ad
import infercnvpy as cnv
import numpy as np
import pandas as pd
import scanpy as sc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'comparators'))
from run_infercnvpy import window_table  # noqa: E402

MIN_CPM = 1.0


def main(inp, out):
    t0 = time.time()
    inp, out = Path(inp), Path(out)
    out.mkdir(parents=True, exist_ok=False)
    X = pd.read_csv(inp / 'counts.tsv', sep='\t', index_col=0)
    s = pd.read_csv(inp / 'samples.tsv', sep='\t').set_index('sample').loc[X.columns]
    pos = pd.read_csv(inp / 'positions.tsv', sep='\t').set_index('gene')
    pos['chromosome'] = 'chr' + pos.chromosome.astype(str)
    cpm = X / X.sum(axis=0) * 1e6
    normal = s.role.eq('normal').to_numpy()
    keep = (cpm.loc[:, normal].mean(axis=1) >= MIN_CPM) & cpm.index.isin(pos.index)
    X = X[keep]
    adata = ad.AnnData(X=X.T.to_numpy(np.float32), obs=s[['role']].copy(), var=pos.loc[X.index, ['chromosome', 'start', 'end']])
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    cnv.tl.infercnv(adata, reference_key='role', reference_cat='normal', window_size=100, step=10,
                    n_jobs=int(os.environ.get('SLURM_CPUS_PER_TASK', '1')))
    Z = adata.obsm['X_cnv']
    Z = Z.toarray() if hasattr(Z, 'toarray') else np.asarray(Z)
    _, genes = window_table(adata.var, adata.uns['cnv']['chr_pos'], Z.shape[1], 100, 10)
    np.savez_compressed(out / 'infercnvpy_windows.npz', Z=Z.astype(np.float32), samples=np.asarray(adata.obs_names, dtype=str))
    genes[['gene', 'bin']].to_csv(out / 'infercnvpy_gene_window.tsv', sep='\t', index=False)
    (out / 'infercnvpy_seconds.txt').write_text(f'{time.time() - t0:.1f}\n')


if __name__ == '__main__':
    main(*sys.argv[1:])
