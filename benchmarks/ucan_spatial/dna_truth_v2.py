"""Schema-2 spatial/DNA supported-interval concordance (U-CAN, EVA PRJEB61514).

DNA is read only from ~/CRC/data/external/eva_cnv/vcf. Expression measures dosage relative to the
tumour's own average, so DNA is expressed as log2(TCN / ploidy) and compared per chromosome arm.
Usage: dna_truth_v2.py HAPLODOSE_OUT SAMPLE OUTDIR [--overlap-fraction FRACTION]

Outputs in OUTDIR:
  arms.tsv        per arm: DNA rel_cn, LOH fraction, min CF; RNA cell-weighted tumour log ratio, dose, state
  subclonal.tsv   continuous callable overlap, chromosome-preserving null and Monte Carlo uncertainty
  summary.json    arm concordance and explicit limits; bulk DNA cannot orient RNA clone contrasts
"""
import gzip
import argparse
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


def union_intervals(intervals):
    out = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if out and start <= out[-1][1]:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([float(start), float(end)])
    return out


def interval_bp(intervals):
    return sum(end - start for start, end in union_intervals(intervals))


def supported_intervals(start, end, support):
    return [(max(start, a), min(end, b)) for a, b in support if min(end, b) > max(start, a)]


def compress_intervals(intervals, support):
    """Map genomic intervals into chromosome-specific callable coordinates."""
    result, offset = [], 0.0
    for a, b in support:
        for start, end in intervals:
            lo, hi = max(start, a), min(end, b)
            if hi > lo:
                result.append((offset + lo - a, offset + hi - a))
        offset += b - a
    return union_intervals(result)


def event_overlap(start, end, calls, support, expected_direction=None):
    """Signed overlap needs independently oriented clone-contrast DNA truth."""
    support = union_intervals(support)
    denominator = interval_bp(supported_intervals(start, end, support))
    hits, signed = [], []
    for row in calls.itertuples():
        overlap = supported_intervals(max(start, row.start), min(end, row.end), support)
        hits.extend(overlap)
        delta = getattr(row, 'mean_diff', np.nan)
        if expected_direction is not None and np.isfinite(delta) and np.sign(delta) == np.sign(expected_direction):
            signed.extend(overlap)
    return dict(callable_event_bp=denominator, overlap_bp=interval_bp(hits),
                overlap_fraction=interval_bp(hits) / denominator if denominator else np.nan,
                direction_evaluable=bool(expected_direction is not None and denominator),
                signed_overlap_fraction=interval_bp(signed) / denominator
                if expected_direction is not None and denominator else np.nan)


def callable_shift_null(start, end, calls, support, rng, draws=200):
    """Preserve chromosome and callable event length, not unassayed physical gaps."""
    support = union_intervals(support)
    size = interval_bp(support)
    length = interval_bp(supported_intervals(start, end, support))
    if size <= 0 or length <= 0:
        return np.full(draws, np.nan)
    cc = compress_intervals(list(zip(calls.start, calls.end)), support)
    values = []
    for origin in rng.uniform(0, size, size=draws):
        end_ = origin + length
        parts = [(origin, min(end_, size))] + ([(0, end_ - size)] if end_ > size else [])
        overlap = [(max(a, c), min(b, d)) for a, b in parts for c, d in cc if min(b, d) > max(a, c)]
        values.append(interval_bp(overlap) / length)
    return np.asarray(values)


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
            rows.append(dict(chromosome=f[0], start=int(f[1]) - 1, end=int(info['END']), logratio=num('CNLR_MEDIAN'),
                             tcn=num('TCN_EM'), lcn=num('LCN_EM'), cf=num('CF_EM')))
    seg = pd.DataFrame(rows)
    return seg, head['purity'], head['ploidy'], vcfs[0]


def arm_of(chrom, pos, cen):
    return chrom + np.where(pos < chrom.map(cen).to_numpy(), 'p', 'q')


def facets_arms(sample, cen):
    """FACETS per arm: length-weighted rel_cn = log2(TCN / ploidy), TCN, LOH fraction, min CF, state."""
    seg, purity, ploidy, vcf = read_facets(sample)
    seg = seg[seg.chromosome.isin(cen.index)].copy()
    seg['len'] = seg.end - seg.start
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


def main(out, sample, dest, overlap_fraction=None):
    if overlap_fraction is not None and not 0 < overlap_fraction <= 1:
        raise ValueError('overlap fraction must be in (0, 1]')
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
    if not len(diff):
        diff = pd.DataFrame(columns=['chromosome', 'start', 'end', 'mean_diff'])
    # RNA bins define observable interval support. Bulk FACETS cannot orient a
    # pair of RNA clones: retain unknown signed truth rather than invent a sign.
    rng = np.random.default_rng(20261101)
    event_rows = []
    for r in sub.itertuples():
        d = diff[diff.chromosome.eq(r.chromosome)]
        support_bins = bins[bins.chromosome.eq(r.chromosome)]
        support = union_intervals(list(zip(support_bins.start, support_bins.end)))
        metric = event_overlap(r.start, r.end, d, support)
        null = callable_shift_null(r.start, r.end, d, support, rng)
        finite = null[np.isfinite(null)]
        metric.update(direction_reason='bulk DNA does not identify oriented RNA clone pairs',
                      signed_truth='unknown',
                      overlap_pass=bool(metric['overlap_fraction'] >= overlap_fraction)
                      if overlap_fraction is not None and metric['callable_event_bp'] else None,
                      null_mean_overlap_fraction=float(finite.mean()) if len(finite) else np.nan,
                      null_mcse=float(finite.std(ddof=1) / np.sqrt(len(finite))) if len(finite) > 1 else np.nan,
                      null_q025=float(np.quantile(finite, 0.025)) if len(finite) else np.nan,
                      null_q975=float(np.quantile(finite, 0.975)) if len(finite) else np.nan)
        event_rows.append(metric)
    if event_rows:
        sub = pd.concat([sub.reset_index(drop=True), pd.DataFrame(event_rows)], axis=1)
    else:
        for col in ['overlap_fraction', 'callable_event_bp', 'direction_evaluable', 'overlap_pass']:
            sub[col] = pd.Series(dtype=float)
    sub.assign(sample=sample).to_csv(dest / 'subclonal.tsv', sep='\t', index=False)

    rho = lambda a, b: float(spearmanr(a, b).statistic) if len(a) > 3 and a.notna().all() else None
    summary = dict(
        schema_version=2, interval_convention='0-based half-open; VCF POS converted from 1-based',
        sample=sample, vcf=str(vcf), purity=purity, ploidy=ploidy, n_arms=int(len(table)),
        tumour_clones=int(len(tum)), tumour_cells=int(tum.sum()),
        rho_log_ratio=rho(table.rna_log_ratio, table.rel_cn), rho_log_dose=rho(table.rna_log_dose, table.rel_cn),
        rho_log_ratio_gene=rho(table.rna_log_ratio_gene, table.rel_cn) if 'rna_log_ratio_gene' in table else None,
        dna_altered_arms=int(altered.sum()),
        altered_arms_same_direction=int((table.rna_state.eq(table.dna_state) & altered).sum()),
        altered_arms_sign_agree=int((np.sign(table.rna_log_ratio) == np.where(table.dna_state.eq('gain'), 1, -1))[altered].sum()),
        neutral_arms=int((~altered).sum()), neutral_arms_called_altered=int((table.rna_state.ne('neutral') & ~altered).sum()),
        subclonal_segments=int(len(sub)),
        subclonal_callable_events=int(sub.callable_event_bp.gt(0).sum()),
        subclonal_mean_overlap_fraction=float(sub.overlap_fraction.mean()) if len(sub) else None,
        direction_evaluable_events=0,
        direction_limit='bulk FACETS cannot orient RNA clone contrasts; unsigned overlap is not clone validation',
        subclonal_overlap_pass=int(sub.overlap_pass.fillna(False).sum()) if overlap_fraction is not None else None,
        null_model='within-chromosome circular shifts in callable coordinates; preserves callable event bp, '
                   'not physical gap lengths; mean overlap and Monte Carlo uncertainty per event',
        thresholds=dict(gain=GAIN, subclonal_cf_over_purity=SUBCLONAL,
                        min_subclonal_bp=MIN_SUBCLONAL_BP, overlap_fraction=overlap_fraction))
    (dest / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Schema-2 supported spatial/DNA interval concordance')
    ap.add_argument('out')
    ap.add_argument('sample')
    ap.add_argument('dest')
    ap.add_argument('--overlap-fraction', type=float, default=None,
                    help='explicit descriptive overlap threshold; default reports continuous overlap only')
    a = ap.parse_args()
    main(a.out, a.sample, a.dest, a.overlap_fraction)
