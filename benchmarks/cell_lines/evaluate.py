"""Score a Haplodose run on the GSE149224 cell-line design (see prepare.py) against DepMap 24Q4 copy number.

DepMap OmicsCNGene is gene-level relative copy number on a log scale; only ranks and signs are used here, so
the exact transform does not matter. Truth for a line is its DepMap value minus RKO's (RKO is the reference).
Usage: evaluate.py HAPLODOSE_OUT INPUT_DIR DEST
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score

DEPMAP = Path('/gpfs/commons/home/jameslee/CRC/data/external/depmap_24q4/OmicsCNGene.csv')
MODELS = {'ACH-000971': 'HCT116', 'ACH-000842': 'SW480', 'ACH-000943': 'RKO'}
CEN = Path(__file__).resolve().parents[2] / 'python/haplodose/data/centromeres_hg38.tsv'
DIFF = 0.2   # |DepMap difference| treated as a real CN difference [scale per DepMap docs, not verified here]


def depmap():
    header = subprocess.run(['head', '-1', str(DEPMAP)], capture_output=True, text=True, check=True).stdout
    rows = subprocess.run(['grep', '-E', '^(' + '|'.join(MODELS) + '),', str(DEPMAP)],
                          capture_output=True, text=True, check=True).stdout.splitlines()
    cols = [c.split(' (')[0] for c in header.rstrip('\n').split(',')[1:]]
    dm = pd.DataFrame([[float(x) if x else np.nan for x in r.split(',')[1:]] for r in rows], columns=cols,
                      index=[MODELS[r.split(',')[0]] for r in rows]).T
    return dm.loc[~dm.index.duplicated()]


def main(out, inp, dest):
    out, dest = Path(out), Path(dest)
    dest.mkdir(parents=True, exist_ok=False)
    dm = depmap()
    truth = pd.DataFrame({l: dm[l] - dm['RKO'] for l in ['HCT116', 'SW480']})
    genes = pd.read_csv(out / 'genes.tsv.gz', sep='\t')
    assert {'gene', 'bin'} <= set(genes.columns), genes.columns
    bins = pd.read_csv(out / 'bins.tsv', sep='\t')
    cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    bins['arm'] = bins.chromosome + np.where(bins.start < bins.chromosome.map(cen), 'p', 'q')
    gt = genes.merge(truth, left_on='gene', right_index=True)
    bin_truth = gt.groupby('bin')[['HCT116', 'SW480']].median()
    cells = pd.read_csv(out / 'cells.tsv.gz', sep='\t')
    prof = pd.read_csv(out / 'profiles.tsv.gz', sep='\t')
    clones = pd.read_csv(out / 'clones.tsv', sep='\t')
    q = cells[cells.design.isin(['mixture', 'single_line'])]
    rows, spec_rows = [], []
    for (spec, line), c in q.groupby(['specimen', 'truth_line']):
        w = c[c.malignant_call].clone.value_counts()
        if not len(w):
            continue
        p = prof[prof.specimen.eq(spec) & prof.clone.isin(w.index)].pivot(index='bin', columns='clone', values='log_ratio')
        r = (p[w.index] * w.to_numpy()).sum(axis=1) / w.sum()
        b = pd.DataFrame({'rna': r}).join(bin_truth[line].rename('dna'), how='inner').join(bins.set_index('bin').arm)
        a = b.groupby('arm')[['rna', 'dna']].median()
        rows.append(dict(specimen=spec, line=line, design=c.design.iloc[0], n_cells=len(c),
                         malignant_fraction=float(c.malignant_call.mean()),
                         rho_bins=float(spearmanr(b.rna, b.dna).statistic), rho_arms=float(spearmanr(a.rna, a.dna).statistic)))
    for spec, c in q.groupby('specimen'):
        tum = clones[clones.specimen.eq(spec) & ~clones.normal_like]
        d = dict(specimen=spec, design=c.design.iloc[0], n_cells=len(c), tumour_clones=len(tum),
                 ari_line=float(adjusted_rand_score(c.truth_line, c.clone)) if c.truth_line.nunique() > 1 else np.nan,
                 ari_dose=float(adjusted_rand_score(c.dose, c.clone)) if c.dose.nunique() > 1 else np.nan,
                 ari_batch=float(adjusted_rand_score(c.batch, c.clone)) if c.batch.nunique() > 1 else np.nan)
        spec_rows.append(d)
    lines = pd.DataFrame(rows)
    specs = pd.DataFrame(spec_rows)
    # Reported differences between tumour clones in mixtures: sign vs DepMap (majority line of each clone).
    diff = pd.read_csv(out / 'differences.tsv', sep='\t')
    maj = q[q.malignant_call].groupby('clone').truth_line.agg(lambda s: s.value_counts().index[0])
    gb = genes.merge(dm[['HCT116', 'SW480']], left_on='gene', right_index=True)
    drows = []
    for r in diff[diff.specimen.str.startswith('mix_')].itertuples():
        la, lb = maj.get(r.clone_a), maj.get(r.clone_b)
        if la is None or lb is None or la == lb:
            continue
        g = gb[(gb.bin >= r.first_bin) & (gb.bin <= r.last_bin)]
        d = float((g[la] - g[lb]).median())             # clone_a minus clone_b, as mean_diff (profile.difference_segments)
        drows.append(dict(specimen=r.specimen, clone_a=r.clone_a, clone_b=r.clone_b, line_a=la, line_b=lb,
                          chromosome=r.chromosome, n_genes=r.n_genes, mean_diff=r.mean_diff, depmap_diff=d,
                          sign_agrees=bool(np.sign(d) == np.sign(r.mean_diff)), depmap_supported=bool(abs(d) >= DIFF)))
    dd = pd.DataFrame(drows)
    ref = pd.read_csv(out / 'reference_check.tsv', sep='\t')
    lines.to_csv(dest / 'line_profiles.tsv', sep='\t', index=False)
    specs.to_csv(dest / 'specimens.tsv', sep='\t', index=False)
    dd.to_csv(dest / 'differences_vs_depmap.tsv', sep='\t', index=False)
    mixs, singles = specs[specs.design.eq('mixture')], specs[specs.design.eq('single_line')]
    summary = dict(
        depmap_genes=int(len(dm)), genes_scored=int(gt.gene.nunique()),
        rho_arms_median={l: float(lines[lines.line.eq(l)].rho_arms.median()) for l in ['HCT116', 'SW480']},
        rho_bins_median={l: float(lines[lines.line.eq(l)].rho_bins.median()) for l in ['HCT116', 'SW480']},
        mixture_ari_line_median=float(mixs.ari_line.median()), mixture_tumour_clones=mixs.tumour_clones.tolist(),
        single_line_tumour_clones=singles.tumour_clones.tolist(), single_line_ari_dose=singles.ari_dose.tolist(), single_line_ari_batch=singles.ari_batch.tolist(),
        differences_between_lines=int(len(dd)),
        differences_sign_agree=int(dd.sign_agrees.sum()) if len(dd) else 0,
        differences_depmap_supported=int((dd.sign_agrees & dd.depmap_supported).sum()) if len(dd) else 0,
        rko_holdout_called_aneuploid=ref.to_dict('records'))
    (dest / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main(*sys.argv[1:4])
