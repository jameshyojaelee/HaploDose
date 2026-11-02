"""Moorman 2025 (HTAN HTA8) scRNA vs MSK-IMPACT tumour DNA, per chromosome arm.

DNA: data_cna_hg19.seg (log2 ratio vs matched normal, not purity adjusted), hg19 coordinates; arms from the hg19
UCSC cytoBand. RNA: the tool's arms.tsv (hg38 bins, hg38 centromeres). Arms are matched by name, so no liftover.
Specimen match: DNA sample's H5AD_SAMPLE_ID == scRNA specimen ('sample'); otherwise same patient ('patient').
Works on Haplodose output and on comparator output in the same layout.
Usage: dna_truth.py TOOL_OUT DEST
"""
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

D = Path('/gpfs/commons/home/jameslee/CRC/data/external/igf2_multimodal_atlas/Moorman_DNA')
CYTO_HG19 = Path('/gpfs/commons/home/jameslee/CRC/data/external/bian_cnv_annotations/cytoBand.txt.gz')
ALTERED = 0.2          # |arm log2 ratio| at or above this is a DNA gain or loss (panel log2, purity diluted)


def hg19_p_arm_end():
    with gzip.open(CYTO_HG19, 'rt') as fh:
        c = pd.read_csv(fh, sep='\t', header=None, names=['chromosome', 'start', 'end', 'band', 'stain'])
    return c[c.stain.eq('acen') & c.band.str.startswith('p')].set_index('chromosome').end


def dna_arms():
    seg = pd.read_csv(D / 'data_cna_hg19.seg', sep='\t')
    seg['chromosome'] = 'chr' + seg.chrom.astype(str).str.replace('^chr', '', regex=True)
    pend = hg19_p_arm_end()
    seg = seg[seg.chromosome.isin(pend.index)]
    parts = []
    for side in ('p', 'q'):
        c = seg.chromosome.map(pend)
        s = seg.assign(start=np.maximum(seg['loc.start'], c) if side == 'q' else seg['loc.start'],
                       end=np.minimum(seg['loc.end'], c) if side == 'p' else seg['loc.end'])
        parts.append(s[s.end > s.start].assign(arm=lambda x: x.chromosome + side))
    s = pd.concat(parts)
    s['w'] = s.end - s.start
    return s.groupby(['ID', 'arm']).apply(lambda g: np.average(g['seg.mean'], weights=g.w),
                                          include_groups=False).rename('dna_log2').reset_index()


def main(out, dest):
    out, dest = Path(out), Path(dest)
    dest.mkdir(parents=True, exist_ok=False)
    samples = pd.read_csv(D / 'data_clinical_sample.txt', sep='\t', comment='#')
    dna = dna_arms().merge(samples[['SAMPLE_ID', 'PATIENT_ID', 'H5AD_SAMPLE_ID', 'SAMPLE_TYPE']],
                           left_on='ID', right_on='SAMPLE_ID')
    dna = dna[~dna.SAMPLE_TYPE.eq('Normal')]
    dna['patient'] = dna.H5AD_SAMPLE_ID.str.extract(r'^(KG\d+)', expand=False)
    clones = pd.read_csv(out / 'clones.tsv', sep='\t')
    arms = pd.read_csv(out / 'arms.tsv', sep='\t')
    prof = pd.read_csv(out / 'profiles.tsv.gz', sep='\t', usecols=['specimen', 'clone', 'bin', 'cn_state'])
    bins = pd.read_csv(out / 'bins.tsv', sep='\t')
    cen = pd.read_csv(Path(__file__).resolve().parents[2] / 'python/haplodose/data/centromeres_hg38.tsv',
                      sep='\t', comment='#').set_index('chromosome').p_arm_end
    bins['arm'] = bins.chromosome + np.where(bins.start < bins.chromosome.map(cen), 'p', 'q')
    rows, arm_rows = [], []
    for spec, tum in clones[~clones.normal_like.astype(bool)].groupby('specimen'):
        w = tum.set_index('clone').n_cells_map.astype(float)
        if w.sum() == 0:
            continue
        a = arms[arms.specimen.eq(spec) & arms.clone.isin(w.index)]
        rna = a.assign(w=a.clone.map(w)).groupby('arm').apply(
            lambda g: np.average(g.median_log_ratio, weights=g.w), include_groups=False).rename('rna_log_ratio')
        p = prof[prof.specimen.eq(spec) & prof.clone.isin(w.index)].merge(bins[['bin', 'arm']], on='bin')
        p['w'] = p.clone.map(w)
        fr = p.assign(g=p.cn_state.ge(3) * p.w, l=p.cn_state.eq(1) * p.w).groupby('arm')[['g', 'l', 'w']].sum()
        state = pd.Series(np.select([fr.g / fr.w > 0.5, fr.l / fr.w > 0.5], ['gain', 'loss'], 'neutral'),
                          index=fr.index, name='rna_state')
        exact = dna[dna.H5AD_SAMPLE_ID.eq(spec)]
        match = 'sample'
        if exact.empty:
            pat = pd.Series([spec]).str.extract(r'^(KG\d+)')[0][0]
            exact = dna[dna.patient.eq(pat)]
            match = 'patient'
            if exact.empty:
                continue
        for sid, d in exact.groupby('SAMPLE_ID'):
            t = d.set_index('arm')[['dna_log2']].join(rna, how='inner').join(state, how='inner')
            t['dna_state'] = np.select([t.dna_log2 >= ALTERED, t.dna_log2 <= -ALTERED], ['gain', 'loss'], 'neutral')
            alt = t.dna_state.ne('neutral')
            arm_rows.append(t.reset_index().assign(specimen=spec, dna_sample=sid, match=match))
            rows.append(dict(specimen=spec, dna_sample=sid, dna_h5ad=d.H5AD_SAMPLE_ID.iloc[0], match=match,
                             tumour_clones=int(len(w)), tumour_cells=int(w.sum()), n_arms=int(len(t)),
                             rho=float(spearmanr(t.rna_log_ratio, t.dna_log2).statistic) if len(t) > 3 else np.nan,
                             dna_altered_arms=int(alt.sum()),
                             altered_same_direction=int((t.rna_state.eq(t.dna_state) & alt).sum()),
                             neutral_arms=int((~alt).sum()),
                             neutral_called_altered=int((t.rna_state.ne('neutral') & ~alt).sum())))
    res = pd.DataFrame(rows)
    res.to_csv(dest / 'specimens.tsv', sep='\t', index=False)
    if arm_rows:
        pd.concat(arm_rows).to_csv(dest / 'arms.tsv', sep='\t', index=False)
    exact = res[res.match.eq('sample')] if len(res) else res
    summary = dict(n_pairs=int(len(res)), n_exact=int(len(exact)),
                   n_exact_patients=int(exact.specimen.str.extract(r'^(KG\d+)')[0].nunique()) if len(exact) else 0,
                   median_rho_exact=float(exact.rho.median()) if len(exact) else None,
                   altered_same_direction=f"{int(exact.altered_same_direction.sum())}/{int(exact.dna_altered_arms.sum())}" if len(exact) else None,
                   neutral_called_altered=f"{int(exact.neutral_called_altered.sum())}/{int(exact.neutral_arms.sum())}" if len(exact) else None,
                   altered_threshold_log2=ALTERED)
    (dest / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(res.to_string(index=False))
    print(json.dumps(summary))


if __name__ == '__main__':
    main(*sys.argv[1:3])
