"""Separate cell selection from profile estimation on the U-CAN Visium HD runs.

For each sample, pool Haplodose's per-cell bin counts (bin_counts.npy) over different cell sets,
take the log ratio to the pooled normal-epithelium reference, and score arm medians against FACETS
(dna/arms.tsv rel_cn, Spearman). If a simple pool over infercnvpy's malignant cells beats Haplodose's
own clone profiles, the loss is in cell selection or clone weighting, not in the bin counts.
Usage: pool_diagnostic.py HAPLODOSE_ROOT INFERCNVPY_ROOT OUT_TSV
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from dna_truth import CEN, arm_of  # noqa: E402


def arm_profile(Y, rows, ref_rows, bin_arm):
    q = Y[rows].sum(axis=0, dtype=np.float64)
    r = Y[ref_rows].sum(axis=0, dtype=np.float64)
    lr = np.log((q + 0.5) / q.sum()) - np.log((r + 0.5) / r.sum())
    lr -= np.median(lr)
    return pd.Series(lr).groupby(bin_arm).median()


def main(hd_root, icp_root, out_tsv):
    hd_root, icp_root = Path(hd_root), Path(icp_root)
    cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    rows = []
    for sdir in sorted(p for p in hd_root.iterdir() if (p / 'dna/arms.tsv').exists()):
        s = sdir.name
        hd = sdir / 'haplodose'
        icp = sorted(icp_root.glob(f'{s}/infercnvpy_*'))
        dna = pd.read_csv(sdir / 'dna/arms.tsv', sep='\t').set_index('arm')
        cells = pd.read_csv(hd / 'cells.tsv.gz', sep='\t',
                            usecols=['cell_id', 'role', 'ref_group', 'malignant_call', 'total_umi_included'])
        bins = pd.read_csv(hd / 'bins.tsv', sep='\t')
        bin_arm = arm_of(bins.chromosome, bins.start.to_numpy(), cen)
        Y = np.load(hd / 'bin_counts.npy', mmap_mode='r')
        query = cells.role.eq('tumour_epithelium').to_numpy()
        ref = (cells.role.eq('reference') & cells.ref_group.fillna('').str.contains('normal_epi')).to_numpy()
        sets = {'all_query': query, 'hd_malignant': query & cells.malignant_call.fillna(False).astype(bool).to_numpy()}
        if icp:
            ic = pd.read_csv(icp[0] / 'cells.tsv', sep='\t').set_index('cell_id').malignant_call
            icp_mal = cells.cell_id.map(ic).fillna(False).astype(bool).to_numpy()
            sets['icp_malignant'] = query & icp_mal
            sets['icp_normal'] = query & ~icp_mal
            sets['hd_mal_and_icp_mal'] = sets['hd_malignant'] & icp_mal
        deep = cells.total_umi_included.to_numpy() >= np.median(cells.total_umi_included[query])
        sets['all_query_deep'] = query & deep
        for name, m in sets.items():
            idx = np.flatnonzero(m)
            if len(idx) < 50:
                continue
            a = arm_profile(Y, idx, np.flatnonzero(ref), bin_arm).reindex(dna.index)
            ok = a.notna()
            b, _ = np.polyfit(dna.rel_cn[ok], a[ok], 1)
            rows.append(dict(sample=s, cells=name, n=len(idx), rho=spearmanr(a[ok], dna.rel_cn[ok])[0], slope=b))
        rows.append(dict(sample=s, cells='hd_reported', n=np.nan,
                         rho=spearmanr(dna.rna_log_ratio, dna.rel_cn, nan_policy='omit')[0]))
        if icp and (icp[0] / 'dna/arms.tsv').exists():
            d2 = pd.read_csv(icp[0] / 'dna/arms.tsv', sep='\t')
            rows.append(dict(sample=s, cells='icp_reported', n=np.nan,
                             rho=spearmanr(d2.rna_log_ratio, d2.rel_cn, nan_policy='omit')[0]))
        print(s, flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(out_tsv, sep='\t', index=False)
    print(res.pivot(index='sample', columns='cells', values='rho').round(2).to_string())
    print(res.groupby('cells')[['rho', 'slope']].median().round(3).to_string())


if __name__ == '__main__':
    main(*sys.argv[1:])
