"""Do in-section reference groups (immune, stromal bins of the same Visium HD section) carry tumour copy number?

Per sample: pool each reference group's bins, take the log ratio to the pooled external normal-epithelium
reference, and correlate arm medians with FACETS rel_cn (Spearman). Also lists the fitted baseline mixture
weights from Haplodose's specimens.tsv. A group whose profile tracks the tumour DNA leaks tumour RNA into the
baseline and, with a high mixture weight, cancels the copy-number signal.
Usage: reference_leak_check.py HAPLODOSE_ROOT OUT_TSV
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from dna_truth import CEN, arm_of  # noqa: E402


def main(root, out_tsv):
    cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    res = []
    for sdir in sorted(p for p in Path(root).iterdir() if (p / 'dna/arms.tsv').exists()):
        hd = sdir / 'haplodose'
        dna = pd.read_csv(sdir / 'dna/arms.tsv', sep='\t').set_index('arm')
        cells = pd.read_csv(hd / 'cells.tsv.gz', sep='\t', usecols=['role', 'ref_group'])
        bins = pd.read_csv(hd / 'bins.tsv', sep='\t')
        arm = arm_of(bins.chromosome, bins.start.to_numpy(), cen)
        Y = np.load(hd / 'bin_counts.npy', mmap_mode='r')
        mix = dict(kv.split('=') for kv in pd.read_csv(hd / 'specimens.tsv', sep='\t').mixture[0].split(';'))
        ref = cells.role.eq('reference').to_numpy()
        epi = ref & cells.ref_group.fillna('').str.contains('normal_epi').to_numpy()
        r = Y[np.flatnonzero(epi)].sum(axis=0, dtype=np.float64)
        groups = {'query': cells.role.eq('tumour_epithelium').to_numpy()}
        groups.update({g: ref & cells.ref_group.eq(g).to_numpy() for g in ('immune', 'stromal')})
        row = dict(sample=sdir.name, **{f'w_{k.split("|")[1]}': float(v) for k, v in mix.items()
                                        if 'normal_epi' not in k},
                   w_normal_epi=sum(float(v) for k, v in mix.items() if 'normal_epi' in k))
        for g, m in groups.items():
            q = Y[np.flatnonzero(m)].sum(axis=0, dtype=np.float64)
            lr = np.log((q + 0.5) / q.sum()) - np.log((r + 0.5) / r.sum())
            a = pd.Series(lr - np.median(lr)).groupby(arm).median().reindex(dna.index)
            ok = a.notna()
            row[f'rho_{g}_vs_dna'] = spearmanr(a[ok], dna.rel_cn[ok])[0]
            row[f'n_{g}'] = int(m.sum())
        row['rho_hd_reported'] = spearmanr(dna.rna_log_ratio, dna.rel_cn, nan_policy='omit')[0]
        res.append(row)
    res = pd.DataFrame(res)
    res.to_csv(out_tsv, sep='\t', index=False)
    print(res.round(2).to_string())


if __name__ == '__main__':
    main(*sys.argv[1:])
