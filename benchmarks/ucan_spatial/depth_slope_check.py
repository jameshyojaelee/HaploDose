"""Does the depth slope absorb copy number on Visium HD?

Re-derives Haplodose's initial per-bin depth slope (SpecimenModel.depth_slopes with one clone: depth
quintiles, log share regressed on mean log depth) from bin_counts.npy, once over the query tiles and
once over the normal-epithelium reference tiles, which carry no tumour copy number. Arm medians of each
slope are compared with FACETS rel_cn (Spearman). A query slope that tracks DNA while the reference
slope does not means the depth term is learning tumour content, not depth.
Usage: depth_slope_check.py HAPLODOSE_ROOT OUT_TSV
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from dna_truth import CEN, arm_of  # noqa: E402


def slopes(Y, rows, n_q=5):
    L = Y[rows].sum(axis=1, dtype=np.float64)
    x = np.log(np.maximum(L, 1.0))
    q = np.searchsorted(np.quantile(x, np.linspace(0, 1, n_q + 1)[1:-1]), x)
    share, xg = [], []
    for g in range(n_q):
        r = rows[q == g]
        y = Y[r].sum(axis=0, dtype=np.float64)
        share.append(np.log((y + 0.5) / y.sum()))
        xg.append(np.average(x[q == g], weights=L[q == g]))
    share, xg = np.array(share), np.array(xg)
    share -= share.mean(axis=0)
    xc = xg - xg.mean()
    return xc @ share / (xc @ xc)


def main(root, out_tsv):
    cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    res = []
    for sdir in sorted(p for p in Path(root).iterdir() if (p / 'dna/arms.tsv').exists()):
        hd = sdir / 'haplodose'
        dna = pd.read_csv(sdir / 'dna/arms.tsv', sep='\t').set_index('arm')
        cells = pd.read_csv(hd / 'cells.tsv.gz', sep='\t', usecols=['role', 'ref_group', 'malignant_call'])
        bins = pd.read_csv(hd / 'bins.tsv', sep='\t')
        arm = arm_of(bins.chromosome, bins.start.to_numpy(), cen)
        Y = np.load(hd / 'bin_counts.npy', mmap_mode='r')
        query = np.flatnonzero(cells.role.eq('tumour_epithelium'))
        ref = np.flatnonzero(cells.role.eq('reference') & cells.ref_group.fillna('').str.contains('normal_epi'))
        row = dict(sample=sdir.name)
        for name, rows in (('query', query), ('reference', ref)):
            g = pd.Series(slopes(Y, rows)).groupby(arm).median().reindex(dna.index)
            ok = g.notna()
            row[f'rho_{name}_slope_vs_dna'] = spearmanr(g[ok], dna.rel_cn[ok])[0]
            row[f'sd_{name}_arm_slope'] = g.std()
        row['rho_hd_reported'] = spearmanr(dna.rna_log_ratio, dna.rel_cn, nan_policy='omit')[0]
        res.append(row)
        print(row, flush=True)
    res = pd.DataFrame(res)
    res.to_csv(out_tsv, sep='\t', index=False)
    print(res.round(2).to_string())
    print(res.drop(columns='sample').median().round(3).to_string())


if __name__ == '__main__':
    main(*sys.argv[1:])
