"""Compare one Haplodose Visium HD run with the patient's WGS FACETS calls (U-CAN, EVA PRJEB61514).

DNA is read only from ~/CRC/data/external/eva_cnv/vcf. Expression measures dosage relative to the
tumour's own average, so DNA is expressed as log2(TCN / ploidy) and compared per chromosome arm.
Usage: dna_truth.py HAPLODOSE_OUT SAMPLE OUTDIR

Outputs in OUTDIR:
  arms.tsv        per arm: DNA rel_cn, LOH fraction, min CF; RNA cell-weighted tumour log ratio, dose, state
  subclonal.tsv   DNA subclonal altered segments (>= 10 Mb) and whether Haplodose split tumour clones there
  summary.json    Spearman rho, arm direction concordance, subclonal recall
"""
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

VCF_DIR = Path('/gpfs/commons/home/jameslee/CRC/data/external/eva_cnv/vcf')
CEN = Path(__file__).resolve().parents[2] / 'python/haplodose/data/centromeres_hg38.tsv'
GAIN = 0.5            # TCN - ploidy at or above this is a DNA gain (below -GAIN a loss)
SUBCLONAL = 0.8       # CF below this fraction of purity is subclonal
MIN_SUBCLONAL_BP = 10_000_000


def read_facets(sample):
    vcfs = list(VCF_DIR.glob(f'*CRC-SW-{sample}-T*.vcf.gz'))
    if len(vcfs) != 1:
        raise FileNotFoundError(f'need one FACETS VCF for {sample}, found {len(vcfs)}')
    head, rows = {}, []
    with gzip.open(vcfs[0], 'rt') as fh:
        for line in fh:
            if line.startswith('##purity=') or line.startswith('##ploidy='):
                k, v = line[2:].strip().split('=')
                head[k] = float(v)
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            info = dict(x.split('=', 1) for x in f[7].split(';') if '=' in x)
            num = lambda k: float(info[k]) if info.get(k, '.') != '.' else np.nan
            rows.append(dict(chromosome=f[0], start=int(f[1]), end=int(info['END']), logratio=num('CNLR_MEDIAN'),
                             tcn=num('TCN_EM'), lcn=num('LCN_EM'), cf=num('CF_EM')))
    seg = pd.DataFrame(rows)
    return seg, head['purity'], head['ploidy'], vcfs[0]


def arm_of(chrom, pos, cen):
    return chrom + np.where(pos < chrom.map(cen).to_numpy(), 'p', 'q')


def facets_arms(sample, cen):
    """FACETS per arm: length-weighted rel_cn = log2(TCN / ploidy), TCN, LOH fraction, min CF, state."""
    seg, purity, ploidy, vcf = read_facets(sample)
    seg = seg[seg.chromosome.isin(cen.index)].copy()
    seg['len'] = seg.end - seg.start + 1
    seg['rel_cn'] = np.log2(np.maximum(seg.tcn, 0.25) / ploidy)
    # Expected minor-allele fraction of the bulk sample (tumour at purity plus 1:1 normal cells).
    seg['minor_expected'] = (purity * seg.lcn + 1 - purity) / (purity * seg.tcn + 2 * (1 - purity))

    # DNA per arm: split segments at the centromere, length-weighted.
    parts = []
    for side in ('p', 'q'):
        c = seg.chromosome.map(cen)
        s = seg.assign(start=np.maximum(seg.start, c) if side == 'q' else seg.start,
                       end=np.minimum(seg.end, c) if side == 'p' else seg.end)
        s = s[s.end > s.start].assign(arm=lambda d: d.chromosome + side)
        parts.append(s)
    s = pd.concat(parts)
    s['w'] = s.end - s.start
    ok = s.tcn.notna()
    dna = s[ok].groupby('arm').apply(lambda g: pd.Series(dict(
        rel_cn=np.average(g.rel_cn, weights=g.w), tcn=np.average(g.tcn, weights=g.w),
        loh_fraction=np.average(g.lcn.eq(0), weights=g.w), dna_bp=g.w.sum(),
        minor_expected=np.average(g.minor_expected[g.lcn.notna()], weights=g.w[g.lcn.notna()])
        if g.lcn.notna().any() else np.nan,
        imbalance_fraction=np.average((g.tcn - 2 * g.lcn).abs().ge(1), weights=g.w) if g.lcn.notna().any() else np.nan,
        min_cf_altered=g.loc[(g.tcn != 2) | (g.lcn != 1), 'cf'].min())), include_groups=False)
    dna['dna_state'] = np.select([dna.tcn - ploidy >= GAIN, dna.tcn - ploidy <= -GAIN], ['gain', 'loss'], 'neutral')
    return dna, seg, purity, ploidy, vcf


def main(out, sample, dest):
    out, dest = Path(out), Path(dest)
    dest.mkdir(parents=True, exist_ok=False)
    cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    dna, seg, purity, ploidy, vcf = facets_arms(sample, cen)

    # RNA per arm: tumour clones only, weighted by MAP cell counts.
    clones = pd.read_csv(out / 'clones.tsv', sep='\t')
    tum = clones[~clones.normal_like].set_index('clone').n_cells_map
    arms = pd.read_csv(out / 'arms.tsv', sep='\t')
    arms = arms[arms.clone.isin(tum.index)]
    w = arms.clone.map(tum).astype(float)
    rna = arms.assign(w=w).groupby('arm').apply(lambda g: pd.Series(dict(
        rna_log_ratio=np.average(g.median_log_ratio, weights=g.w),
        rna_log_dose=np.average(g.mean_log_dose, weights=g.w), rna_bins=g.n_bins.iloc[0],
        # Gene-level column (added 2026-11-02); NaN for runs that predate it.
        rna_log_ratio_gene=np.average(g.median_log_ratio_gene, weights=g.w) if 'median_log_ratio_gene' in g
        else np.nan)), include_groups=False)
    prof = pd.read_csv(out / 'profiles.tsv.gz', sep='\t', usecols=['clone', 'bin', 'cn_state'])
    bins = pd.read_csv(out / 'bins.tsv', sep='\t')
    bins['arm'] = arm_of(bins.chromosome, bins.start.to_numpy(), cen)
    prof = prof[prof.clone.isin(tum.index)].merge(bins[['bin', 'arm']], on='bin')
    prof['w'] = prof.clone.map(tum).astype(float)
    # cn_state 1 loss, 2 neutral, 3 gain, 4 amplification (pipeline.summarise: path + 1)
    frac = prof.assign(gain=prof.cn_state.ge(3) * prof.w, loss=prof.cn_state.eq(1) * prof.w).groupby('arm')[
        ['gain', 'loss', 'w']].sum()
    rna['rna_gain_frac'] = frac.gain / frac.w
    rna['rna_loss_frac'] = frac.loss / frac.w
    rna['rna_state'] = np.select([rna.rna_gain_frac > 0.5, rna.rna_loss_frac > 0.5], ['gain', 'loss'], 'neutral')

    if tum.empty:      # no tumour clone found: every arm is a miss; rho undefined
        rna = pd.DataFrame(index=dna.index, data=dict(rna_log_ratio=np.nan, rna_log_dose=np.nan, rna_bins=np.nan,
                                                      rna_gain_frac=0.0, rna_loss_frac=0.0, rna_state='neutral'))
    table = dna.join(rna, how='inner').reset_index().assign(sample=sample)
    table.to_csv(dest / 'arms.tsv', sep='\t', index=False)
    altered = table.dna_state.ne('neutral')

    # Subclonal DNA events and whether Haplodose split tumour clones on them.
    sub = seg[seg.tcn.notna() & ((seg.tcn != 2) | (seg.lcn != 1)) & (seg.cf < SUBCLONAL * purity)
              & (seg.len >= MIN_SUBCLONAL_BP)].copy()
    try:
        diff = pd.read_csv(out / 'differences.tsv', sep='\t')
    except pd.errors.EmptyDataError:
        diff = pd.DataFrame()
    if len(diff):
        diff = diff[diff.clone_a.isin(tum.index) & diff.clone_b.isin(tum.index)]
    def hit(r):
        if not len(diff):
            return False
        d = diff[diff.chromosome.eq(r.chromosome)]
        return bool(((d.start <= r.end) & (d.end >= r.start)).any())
    sub['tumour_clone_split'] = [hit(r) for r in sub.itertuples()] if len(sub) else []
    # Chance level: the same segments placed at random positions on random chromosomes (length kept,
    # chromosome length from the FACETS segments), 200 draws; mean fraction that hits a reported difference.
    rng = np.random.default_rng(20261101)
    clen = seg.groupby('chromosome').end.max()
    clen = clen[clen.index.isin(cen.index)]
    expected = np.nan
    if len(sub):
        hits = []
        for _ in range(200):
            k = 0
            for r in sub.itertuples():
                ok = clen[clen > r.len]
                c = rng.choice(ok.index)
                st = int(rng.integers(0, ok[c] - r.len))
                k += hit(pd.Series(dict(chromosome=c, start=st, end=st + r.len)))
            hits.append(k / len(sub))
        expected = float(np.mean(hits))
    sub.assign(sample=sample).to_csv(dest / 'subclonal.tsv', sep='\t', index=False)

    rho = lambda a, b: float(spearmanr(a, b).statistic) if len(a) > 3 and a.notna().all() else None
    summary = dict(
        sample=sample, vcf=str(vcf), purity=purity, ploidy=ploidy, n_arms=int(len(table)),
        tumour_clones=int(len(tum)), tumour_cells=int(tum.sum()),
        rho_log_ratio=rho(table.rna_log_ratio, table.rel_cn), rho_log_dose=rho(table.rna_log_dose, table.rel_cn),
        rho_log_ratio_gene=rho(table.rna_log_ratio_gene, table.rel_cn) if 'rna_log_ratio_gene' in table else None,
        dna_altered_arms=int(altered.sum()),
        altered_arms_same_direction=int((table.rna_state.eq(table.dna_state) & altered).sum()),
        altered_arms_sign_agree=int((np.sign(table.rna_log_ratio) == np.where(table.dna_state.eq('gain'), 1, -1))[altered].sum()),
        neutral_arms=int((~altered).sum()), neutral_arms_called_altered=int((table.rna_state.ne('neutral') & ~altered).sum()),
        subclonal_segments=int(len(sub)), subclonal_segments_split=int(sub.tumour_clone_split.sum()) if len(sub) else 0,
        subclonal_split_fraction_expected_by_chance=expected,
        thresholds=dict(gain=GAIN, subclonal_cf_over_purity=SUBCLONAL, min_subclonal_bp=MIN_SUBCLONAL_BP))
    (dest / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    main(*sys.argv[1:4])
