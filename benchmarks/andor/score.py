#!/usr/bin/env python3
"""Score one tool run on one Andor 2020 line (scRNA, Haplodose layout) against that line's scDNA clones.

Usage: score.py TOOL_OUT SCDNA_DIR LINE DEST [ARM_COLUMN]
ARM_COLUMN: arms.tsv column for the RNA arm values (default median_log_ratio; Haplodose also writes
median_log_ratio_gene).
SCDNA_DIR is results/andor_<stamp>/scdna (scdna_cn.py output, one subdirectory per line).
Scored per query specimen (scRNA library) of LINE. RNA clones = clones of the specimen's query cells
with >= max(20, 2% of cells) cells; normal_like clones are reported but kept (a cell line has no
normal cells, so a normal_like clone means "close to the other-line reference").

1. Clone counts: RNA clones vs scDNA clones (scdna qc.json).
2. Clone-difference test (reference-free: the other-line reference cancels between clones of the
   same line). Truth per arm: 1 if any pair of scDNA clones differs on >= 30% of the arm's bins by
   >= 1 copy (clone_differences.tsv). RNA score per arm: largest |difference in median_log_ratio|
   over pairs of RNA clones. Reported: AUROC of the RNA score for the truth (0.5 = chance), and
   precision/recall of the arms hit by the tool's own differences.tsv calls (Haplodose: statistical
   test; comparators: cn_state runs), with the truth base rate as the chance precision.
3. Line profile, matched contrast. The RNA arm profile is relative to the other lines (the
   reference), so DNA is put on the same contrast: per arm, log2(CN_line / ploidy_line) minus the
   mean over the other lines of log2(CN / ploidy), each line's CN the cell-weighted mean of its
   scDNA clones and ploidy its genome-wide mean CN. Spearman over arms. This is not an absolute
   profile check.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

CEN = Path(__file__).resolve().parents[2] / 'python/haplodose/data/centromeres_hg38.tsv'


def arm_of(chrom, start):
    cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    chrom = pd.Series(np.asarray(chrom, dtype=str))
    return (chrom + np.where(np.asarray(start, float) >= chrom.map(cen).fillna(np.inf).to_numpy(), 'q', 'p')).to_numpy()


def dna_line(scdna, line):
    d = Path(scdna) / line
    qc = json.loads((d / 'qc.json').read_text())
    arms = pd.read_csv(d / 'clone_arms.tsv', sep='\t').set_index('arm')
    w = pd.Series(qc['clones'], dtype=float)
    cons = (arms[w.index] * w).sum(axis=1) / w.sum()
    bins = pd.read_csv(d / 'clone_bins.tsv', sep='\t')
    ploidy = float((bins[w.index] * w).sum(axis=1).mean() / w.sum())
    diffs = pd.read_csv(d / 'clone_differences.tsv', sep='\t')
    return qc, cons, ploidy, set(diffs.arm)


def main(tool_out, scdna, line, dest, arm_column='median_log_ratio'):
    tool_out, dest = Path(tool_out), Path(dest)
    if dest.exists() and any(dest.iterdir()):
        raise FileExistsError(dest)
    dest.mkdir(parents=True, exist_ok=True)
    qc, cons, ploidy, diff_arms = dna_line(scdna, line)
    others = [p.name for p in Path(scdna).iterdir() if p.is_dir() and p.name != line and (p / 'qc.json').exists()]
    if others:
        other = pd.concat([np.log2(np.maximum(dna_line(scdna, o)[1], 0.1) / dna_line(scdna, o)[2]).rename(o) for o in others], axis=1)
        dna_contrast = np.log2(np.maximum(cons, 0.1) / ploidy) - other.mean(axis=1)
    else:   # no other line in SCDNA_DIR: the matched-contrast profile check is NA
        dna_contrast = pd.Series(np.nan, index=cons.index)
    all_arms = cons.index

    cells = pd.read_csv(tool_out / 'cells.tsv.gz', sep='\t', dtype={'clone': str})
    clones = pd.read_csv(tool_out / 'clones.tsv', sep='\t')
    arms = pd.read_csv(tool_out / 'arms.tsv', sep='\t')
    # Haplodose writes a lone newline when no clone pair differs; treat that as no calls.
    diffs = pd.read_csv(tool_out / 'differences.tsv', sep='\t') if (tool_out / 'differences.tsv').read_text().strip() else pd.DataFrame()
    q = cells[cells.role.eq('query') & cells.patient.eq(line)]
    rows = []
    for spec, c in q.groupby('specimen'):
        sizes = c.clone[c.clone.notna() & (c.clone != 'NA')].value_counts()
        keep = sizes[sizes >= max(20, 0.02 * len(c))]
        normal_like = set(clones.loc[clones.normal_like & clones.clone.isin(keep.index), 'clone'])
        a = arms[arms.specimen.eq(spec) & arms.clone.isin(keep.index)].pivot(index='arm', columns='clone', values=arm_column)
        a = a.reindex(all_arms)
        names = list(keep.index)
        score = pd.Series(0.0, index=all_arms)
        for i, x in enumerate(names):
            for y in names[i + 1:]:
                score = np.maximum(score, (a[x] - a[y]).abs().fillna(0))
        truth = pd.Series(all_arms.isin(diff_arms).astype(int), index=all_arms)
        auroc = float(roc_auc_score(truth, score)) if 0 < truth.sum() < len(truth) and len(names) > 1 else None
        called = set()
        if len(diffs):
            dd = diffs[diffs.specimen.eq(spec) & diffs.clone_a.isin(names) & diffs.clone_b.isin(names)]
            for r in dd.itertuples():
                called |= set(arm_of([r.chromosome, r.chromosome], [r.start, r.end]))
        tp = len(called & diff_arms)
        prof = (a[names] * keep[names]).sum(axis=1) / keep[names].sum() if names else pd.Series(np.nan, index=all_arms)
        j = pd.concat([prof.rename('rna'), dna_contrast.rename('dna')], axis=1).dropna()
        rows.append(dict(line=line, specimen=spec, arm_column=arm_column, n_cells=len(c), rna_clones=len(names), rna_normal_like_clones=len(normal_like),
                         dna_clones=len(qc['clones']), dna_cells=qc['n_pass_mapd'],
                         dna_subclonal_arms=len(diff_arms), arms=len(all_arms),
                         diff_auroc=auroc, diff_called_arms=len(called), diff_called_true=tp,
                         diff_precision=tp / len(called) if called else None,
                         diff_recall=tp / len(diff_arms) if diff_arms else None,
                         diff_chance_precision=len(diff_arms) / len(all_arms),
                         profile_rho_matched_contrast=float(spearmanr(j.rna, j.dna).statistic) if len(j) > 3 else None,
                         profile_arms=len(j)))
    t = pd.DataFrame(rows)
    t.to_csv(dest / 'scores.tsv', sep='\t', index=False)
    print(t.to_string())


if __name__ == '__main__':
    main(*sys.argv[1:])
