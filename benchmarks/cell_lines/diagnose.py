"""Where is the SW480 signal lost? Arm-level Spearman vs DepMap (line - RKO) for:
  raw     pooled counts per bin, line cells vs RKO fit cells, log ratio minus median (no model)
  model   Haplodose profile log_ratio of the clones the line's cells are assigned to
Usage: diagnose.py HAPLODOSE_OUT"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from evaluate import CEN, depmap  # noqa: E402

out = Path(sys.argv[1])
dm = depmap()
genes = pd.read_csv(out / 'genes.tsv.gz', sep='\t')
bins = pd.read_csv(out / 'bins.tsv', sep='\t')
cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
bins['arm'] = bins.chromosome + np.where(bins.start < bins.chromosome.map(cen), 'p', 'q')
Y = np.load(out / 'bin_counts.npy', mmap_mode='r')
cells = pd.read_csv(out / 'cells.tsv.gz', sep='\t')
assert len(cells) == Y.shape[0], (len(cells), Y.shape)
prof = pd.read_csv(out / 'profiles.tsv.gz', sep='\t')
ref = cells.ref_split.eq('fit').to_numpy() & cells.truth_line.eq('RKO').to_numpy()
yr = np.asarray(Y[ref]).sum(0) + 0.5
for line in ['HCT116', 'SW480']:
    t = (dm[line] - dm['RKO']).rename('dna')
    bt = genes.merge(t, left_on='gene', right_index=True).groupby('bin').dna.median()
    sel = cells.truth_line.eq(line).to_numpy()
    yl = np.asarray(Y[sel]).sum(0) + 0.5
    raw = np.log(yl / yl.sum()) - np.log(yr / yr.sum())
    raw = raw - np.median(raw)
    c = cells[sel & cells.malignant_call.to_numpy()]
    w = c.groupby(['specimen', 'clone']).size()
    p = prof.set_index(['specimen', 'clone', 'bin']).log_ratio
    model = sum(p.loc[s, k].reindex(bins.bin).to_numpy() * n for (s, k), n in w.items()) / w.sum()
    f = pd.DataFrame({'raw': raw, 'model': model, 'arm': bins.arm}, index=bins.bin).join(bt, how='inner')
    a = f.groupby('arm')[['raw', 'model', 'dna']].median()
    print(f"{line}: arms {len(a)} | rho raw {spearmanr(a.raw, a.dna).statistic:.2f} model "
          f"{spearmanr(a.model, a.dna).statistic:.2f} | bins raw {spearmanr(f.raw, f.dna).statistic:.2f} model "
          f"{spearmanr(f.model, f.dna).statistic:.2f} | raw vs model arms {spearmanr(a.raw, a.model).statistic:.2f}")
