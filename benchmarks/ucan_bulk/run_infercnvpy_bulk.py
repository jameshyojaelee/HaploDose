"""infercnvpy 0.6.1 on bulk RNA samples (samples as cells), the usual way inferCNV is applied to bulk.

normalize_total(1e4) + log1p, reference = normal samples, window 100, step 10 (package defaults). Per arm
value = mean of the windows whose centre lies on that arm. Same gene set as naive.py.
Usage: run_infercnvpy_bulk.py INPUT_DIR OUT_ARMS_TSV  (scrna_cnv env)
"""
import os
import sys
from pathlib import Path

import anndata as ad
import infercnvpy as cnv
import numpy as np
import pandas as pd
import scanpy as sc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'comparators'))
sys.path.insert(0, str(HERE.parent / 'ucan_spatial'))
from run_infercnvpy import window_table  # noqa: E402
from dna_truth import CEN, arm_of  # noqa: E402

MIN_CPM = 1.0


def main(inp, out_tsv):
    inp = Path(inp)
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
    bins, _ = window_table(adata.var, adata.uns['cnv']['chr_pos'], Z.shape[1], 100, 10)
    cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    mid = ((bins.start + bins.end) // 2).to_numpy()
    arm = arm_of(bins.chromosome, mid, cen)
    A = pd.DataFrame(Z, index=adata.obs_names).T.groupby(arm).mean().T
    A.index.name = 'sample'
    A.reset_index().melt(id_vars='sample', var_name='arm', value_name='value').to_csv(out_tsv, sep='\t', index=False)


if __name__ == '__main__':
    main(*sys.argv[1:])
