"""Malignant calls and clones from per-cell CNV profiles (infercnvpy, R inferCNV). Needs scanpy.

Neither tool fixes these rules; ours follow the infercnvpy tutorial (Leiden on CNV PCA, cnv_score
per cluster) with a held-out-reference null instead of choosing tumour clusters by eye:
- score per cell = mean squared profile value (Tirosh 2016).
- rule 'cluster' (default): per specimen, Leiden (resolution 1.0) on PCA of the analysed cells
  (query + held-out reference). Cluster statistic = mean over units of the squared cluster-mean
  profile (signal shared by the cluster's cells; per-cell noise averages out). A cluster is
  malignant if its statistic exceeds the largest statistic among clusters from the same Leiden
  procedure run on the held-out reference cells alone (all specimens). The mean per-cell score was
  dropped: it tracks depth noise (smoke check, SMC06-T: 0/189 tumour cells called, 7.6% of
  held-out cells called). Random-subset nulls were too permissive: Leiden groups similar cells, so its
  cluster means vary more than random subsets (smoke 2, U0247: 37% of held-out cells called).
  Held-out reference cells in a malignant cluster are called malignant too; their rate is the
  false-malignant rate in reference_check.tsv.
- rule 'cell': analysed cell malignant if its score exceeds the 99th percentile of held-out scores.
- rule 'all': every query cell is malignant (for designs where the reference is itself tumour, e.g. cell
  lines against other cell lines, so a held-out-reference null says nothing); clones as below.
- clones: per specimen, Leiden (resolution 0.5) on PCA of the malignant cells; non-malignant
  analysed cells form <specimen>_normal.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _leiden(Z, resolution, seed):
    import anndata as ad
    import scanpy as sc
    sub = ad.AnnData(X=np.asarray(Z, dtype=np.float32))
    sc.pp.pca(sub, n_comps=max(2, min(30, Z.shape[0] - 1, Z.shape[1] - 1)), random_state=seed)
    sc.pp.neighbors(sub, n_neighbors=min(15, Z.shape[0] - 1), random_state=seed)
    sc.tl.leiden(sub, resolution=resolution, flavor='igraph', n_iterations=2, random_state=seed)
    return sub.obs.leiden.astype(int).to_numpy()


def call(Z, cells, query, hold, rule='cluster', quantile=0.99, call_resolution=1.0,
         clone_resolution=0.5, min_cells=30, seed=20261101, n_null=1000):
    score = (np.asarray(Z) ** 2).mean(axis=1)
    hold_scores = score[hold]
    rng = np.random.default_rng(seed)
    malignant = np.zeros(len(cells), dtype=bool)
    analysed = query | hold
    spec = cells.specimen.to_numpy()
    cluster_rows = []
    if rule == 'cell':
        malignant = analysed & (score > np.quantile(hold_scores, quantile))
    elif rule == 'all':
        malignant = np.asarray(query, dtype=bool).copy()
    else:
        hrows = np.flatnonzero(hold)
        hlab = _leiden(Z[hrows], call_resolution, seed)
        stat = lambda r: float((np.asarray(Z[r]).mean(axis=0) ** 2).mean())
        cut = max(stat(hrows[hlab == k]) for k in np.unique(hlab))
        for s in pd.unique(spec[analysed]):
            rows = np.flatnonzero((spec == s) & analysed)
            if len(rows) < min_cells:
                continue
            lab = _leiden(Z[rows], call_resolution, seed)
            for k in np.unique(lab):
                m = rows[lab == k]
                cs = stat(m)
                is_mal = cs > cut
                malignant[m] = is_mal
                cluster_rows.append(dict(specimen=s, cluster=int(k), n_cells=len(m), n_holdout=int(hold[m].sum()),
                                         mean_score=float(score[m].mean()), cluster_stat=cs, null_cut=cut, malignant=bool(is_mal)))
    clone = np.full(len(cells), 'NA', dtype=object)
    for s in pd.unique(spec[analysed]):
        in_s = spec == s
        clone[in_s & analysed & ~malignant] = f'{s}_normal'
        rows = np.flatnonzero(in_s & malignant)
        if len(rows) < min_cells:
            clone[rows] = f'{s}_c1'
            continue
        clone[rows] = [f'{s}_c{k + 1}' for k in _leiden(Z[rows], clone_resolution, seed)]
    return score, malignant, clone, pd.DataFrame(cluster_rows)
