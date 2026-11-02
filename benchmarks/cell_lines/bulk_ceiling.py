"""Upper reference for expression-based CN recovery with a mismatched reference: DepMap bulk RNA.

For HCT116 and SW480, compare arm medians of the bulk log expression ratio (line vs RKO, DepMap
OmicsExpressionProteinCodingGenesTPMLogp1) with arm medians of the DepMap CN difference (line - RKO).
Spearman across arms. Same genes and arm definition as evaluate.py (genes in the Haplodose run).
Usage: bulk_ceiling.py HAPLODOSE_OUT DEST_JSON
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from evaluate import CEN, MODELS, depmap  # noqa: E402

EXPR = Path('/gpfs/commons/home/jameslee/CRC/data/external/depmap_24q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv')


def expression():
    header = subprocess.run(['head', '-1', str(EXPR)], capture_output=True, text=True, check=True).stdout
    rows = subprocess.run(['grep', '-E', '^(' + '|'.join(MODELS) + '),', str(EXPR)],
                          capture_output=True, text=True, check=True).stdout.splitlines()
    cols = [c.split(' (')[0] for c in header.rstrip('\n').split(',')[1:]]
    ex = pd.DataFrame([[float(x) if x else np.nan for x in r.split(',')[1:]] for r in rows], columns=cols,
                      index=[MODELS[r.split(',')[0]] for r in rows]).T
    return ex.loc[~ex.index.duplicated()]


def main(out, dest):
    out, dest = Path(out), Path(dest)
    cn, ex = depmap(), expression()
    genes = pd.read_csv(out / 'genes.tsv.gz', sep='\t')
    bins = pd.read_csv(out / 'bins.tsv', sep='\t')
    cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    bins['arm'] = bins.chromosome + np.where(bins.start < bins.chromosome.map(cen), 'p', 'q')
    g = genes.merge(bins[['bin', 'arm']], on='bin').set_index('gene')
    res = {}
    for line in ['HCT116', 'SW480']:
        d = pd.DataFrame({'rna': ex[line] - ex['RKO'], 'dna': cn[line] - cn['RKO']}).join(g.arm, how='inner').dropna()
        d = d[(ex.loc[d.index, ['RKO', line]] > 1).all(axis=1)]          # expressed in both (log TPM+1 > 1)
        a = d.groupby('arm')[['rna', 'dna']].median()
        res[line] = dict(n_genes=int(len(d)), rho_genes=float(spearmanr(d.rna, d.dna).statistic),
                         rho_arms=float(spearmanr(a.rna, a.dna).statistic))
    dest.write_text(json.dumps(res, indent=2) + '\n')
    print(json.dumps(res))


if __name__ == '__main__':
    main(*sys.argv[1:3])
