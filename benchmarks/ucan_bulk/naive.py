"""Naive bulk baseline: per gene log2(CPM + 1) minus the mean over normal samples, centred per sample by its
median, then the median over each arm's genes. Same gene set as Haplodose bulk (protein coding, mean normal
CPM >= MIN_CPM). This is what most bulk RNA copy-number analyses do before any model.
Usage: naive.py INPUT_DIR OUT_ARMS_TSV
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'ucan_spatial'))
from dna_truth import CEN, arm_of  # noqa: E402

MIN_CPM = 1.0


def main(inp, out_tsv):
    inp = Path(inp)
    X = pd.read_csv(inp / 'counts.tsv', sep='\t', index_col=0)
    s = pd.read_csv(inp / 'samples.tsv', sep='\t')
    pos = pd.read_csv(inp / 'positions.tsv', sep='\t').set_index('gene')
    pos['chromosome'] = 'chr' + pos.chromosome.astype(str)
    cpm = X / X.sum(axis=0) * 1e6
    normal = s.set_index('sample').role.loc[X.columns].eq('normal').to_numpy()
    keep = (cpm.loc[:, normal].mean(axis=1) >= MIN_CPM) & cpm.index.isin(pos.index)
    L = np.log2(cpm[keep] + 1)
    R = L.sub(L.loc[:, normal].mean(axis=1), axis=0)
    R = R - R.median(axis=0)
    cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    p = pos.loc[R.index]
    arm = arm_of(p.chromosome, p.start.to_numpy(), cen)
    A = R.groupby(arm).median().T
    A.index.name = 'sample'
    out = A.reset_index().melt(id_vars='sample', var_name='arm', value_name='value')
    out.to_csv(out_tsv, sep='\t', index=False)


if __name__ == '__main__':
    main(*sys.argv[1:])
