"""Which part of the naive bulk statistic gives its better across-patient arm AUC?

Variants of the gene-level log ratio, each centred per sample by its median and summarised by the median over
an arm's genes (same gene set as naive.py):
  naive          log2(CPM + 1) - mean over normals of log2(CPM + 1)          (mean of logs, 1 CPM pseudocount)
  logmean_cpm1   log2(CPM + 1) - log2(mean normal CPM + 1)                  (log of mean, 1 CPM pseudocount)
  meanlog_cnt05  log2(count + 0.5) - log2(L) - mean normal log2((count + 0.5)/L)   (mean of logs, 0.5 count)
  logmean_cnt05  log2((count + 0.5) / (L * mu_g + 0.5)), mu_g pooled normals (Haplodose bulk, gene level)
Usage: diag_reference.py INPUT_DIR OUTDIR   (writes <variant>_arms.tsv)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'ucan_spatial'))
from dna_truth import CEN, arm_of  # noqa: E402


def main(inp, outdir):
    inp, out = Path(inp), Path(outdir)
    out.mkdir(parents=True, exist_ok=False)
    X = pd.read_csv(inp / 'counts.tsv', sep='\t', index_col=0)
    s = pd.read_csv(inp / 'samples.tsv', sep='\t').set_index('sample').loc[X.columns]
    pos = pd.read_csv(inp / 'positions.tsv', sep='\t').set_index('gene')
    pos['chromosome'] = 'chr' + pos.chromosome.astype(str)
    C = X.to_numpy(np.float64)
    L = C.sum(axis=0)
    cpm = C / L * 1e6
    normal = s.role.eq('normal').to_numpy()
    keep = (cpm[:, normal].mean(axis=1) >= 1.0) & X.index.isin(pos.index)
    C, cpm, genes = C[keep], cpm[keep], X.index[keep]
    mu = C[:, normal].sum(axis=1) / L[normal].sum()
    v = {
        'naive': np.log2(cpm + 1) - np.log2(cpm[:, normal] + 1).mean(axis=1, keepdims=True),
        'logmean_cpm1': np.log2(cpm + 1) - np.log2(cpm[:, normal].mean(axis=1, keepdims=True) + 1),
        'meanlog_cnt05': np.log2((C + 0.5) / L) - np.log2((C[:, normal] + 0.5) / L[normal]).mean(axis=1, keepdims=True),
        'logmean_cnt05': np.log2((C + 0.5) / (L * mu[:, None] + 0.5)),
    }
    cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    p = pos.loc[genes]
    arm = arm_of(p.chromosome, p.start.to_numpy(), cen)
    for name, R in v.items():
        R = R - np.median(R, axis=0, keepdims=True)
        A = pd.DataFrame(R, index=genes, columns=X.columns).groupby(arm).median().T
        A.index.name = 'sample'
        A.reset_index().melt(id_vars='sample', var_name='arm', value_name='value').to_csv(
            out / f'{name}_arms.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main(*sys.argv[1:])
