"""CCLE bulk sub-arm benchmark (prespecified SOTA claim 3, docs/STATUS.md): runs benchmarks/tcga_bulk/subarm.py
unchanged except for its DNA inputs, then the prespecified primary test.

Historical version-1 definitions are retained below for provenance. New runs use
version-2 corrected scoring and paired-mean inference, emit retrospective
diagnostics, and never reissue or redefine the original frozen claim.

DNA inputs swapped in (subarm.py otherwise identical: 10 Mb gene-anchored windows, focal = segment shorter than
25% of its arm with |TCN - ploidy| >= 1, DNA-neutral = |TCN - ploidy| < 0.5, thresholds at the 5th/95th
percentile of DNA-neutral windows per group):
  gene positions  hg38 (Ensembl 112, the positions.tsv table); DepMap 24Q4 segments are hg38, so no liftover
  arms           hg38 UCSC cytoBand (md5 6677dd50caf19b0c3626b99ce67ca073, same source as centromeres_hg38.tsv)
  segments       DepMap PureCN SegmentAbsoluteCN; ploidy = length-weighted mean total CN over chr1-22
Tie-break: random (seeded jitter on hd_dose and infercnvpy), so every method sits at the same realised FPR.

Primary (fixed before any run): focal-loss recall at 5% FPR per line, hd_dose vs naive and hd_dose vs infercnvpy,
paired Wilcoxon signed-rank over lines with >= 1 focal-loss window, Holm over these 2 comparisons. Claim holds if
both Holm p < 0.05 and the mean paired difference is > 0. Bootstrap 95% CI of the mean difference (10,000
resamples of lines, seed 20261103), descriptive. Sensitivity (descriptive): the same excluding lines whose DepMap
CN was used as truth in development (dev_lines.tsv).
Amendment 2026-11-02, before the score stage ran (no CCLE result seen): the same primary test within each RNA
library protocol stratum (OmicsProfiles.Stranded of the line's default RNA profile; GTEx is unstranded; 4 of 309
development lines are stranded vs about half of the others), descriptive.
Usage: score.py OUT INPUT_ROOT ICNV_ROOT
"""
import gzip
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

H = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(H / 'benchmarks'))
from hardening.statistics import SCORE_VERSION, adjust_family, paired_mean_comparison
CYTO38 = H / 'data/external/ucsc_hg38/cytoBand.txt.gz'
DM = H / 'data/external/depmap_24q4_ccle'
COMPARATORS = ['naive', 'infercnvpy']
B, SEED = 10_000, 20261103

_spec = importlib.util.spec_from_file_location('subarm', H / 'benchmarks/tcga_bulk/subarm.py')
subarm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(subarm)


def hg38_arms():
    with gzip.open(CYTO38, 'rt') as fh:
        c = pd.read_csv(fh, sep='\t', header=None, names=['chromosome', 'start', 'end', 'band', 'stain'])
    pend = c[c.stain.eq('acen') & c.band.str.startswith('p')].set_index('chromosome').end
    return pend, c.groupby('chromosome').end.max()


def primary(per, dev):
    """Version 2 retrospective scoring; historical claim rules remain in frozen outputs.

    Mean effects, intervals and tests now use the same statistic. Oracle FPR
    thresholds remain descriptive and do not generate a new confirmatory claim.
    """
    out = {}
    for name, p in (('all', per), ('excluding_dev_lines', per[~per['sample'].isin(dev)])):
        rows = []
        for c in COMPARATORS:
            record = paired_mean_comparison(p['recall_loss_hd_dose'], p[f'recall_loss_{c}'],
                                            seed=SEED, n_boot=B, family_size=len(COMPARATORS))
            rows.append(dict(comparator=c, n_lines=record['n'], mean_hd_dose=record['mean_a'],
                             mean_comparator=record['mean_b'], hd_higher=record['a_higher'],
                             comparator_higher=record['b_higher'], **record))
        rows = adjust_family(rows)
        out[name] = dict(score_version=SCORE_VERSION, comparisons=rows,
                         descriptive_rule_met=all(r['superior'] for r in rows), claim_holds=False,
                         evidence_role='retrospective_regression',
                         limitation='same-cohort DNA-oracle thresholds; CI conditions on reference and threshold')
    return out


def main(out, inp, icnv):
    out, inp, icnv = Path(out), Path(inp), Path(icnv)
    seg = pd.read_csv(inp / 'absolute_segments.tsv.gz', sep='\t')
    tab = pd.read_csv(inp / 'absolute_ploidy.tsv', sep='\t', index_col=0)
    pos = subarm_positions(inp)
    subarm.hg19_genes = lambda: pos
    subarm.hg19_arms = hg38_arms
    subarm.absolute = lambda: (tab, seg)
    subarm.GTF37 = 'hg38 positions (Ensembl 112), no liftover; DepMap 24Q4 PureCN absolute segments'
    groups = sorted(d.name for d in inp.iterdir() if (d / 'counts.tsv').exists())
    sys.argv = ['subarm.py', str(out), '--tiebreak', 'random', '--input-scale', 'raw_counts'] + [f'{t}={inp / t}:{icnv / t}' for t in groups]
    subarm.main()
    per = pd.read_csv(out / 'per_patient.tsv.gz', sep='\t')
    dev = set(pd.read_csv(inp / 'dev_lines.tsv', sep='\t').ModelID)
    res = primary(per, dev)
    dp = pd.read_csv(DM / 'OmicsDefaultModelProfiles.csv')
    prof = pd.read_csv(DM / 'OmicsProfiles.csv').set_index('ProfileID').Stranded
    rna = dp[dp.ProfileType.eq('rna')].drop_duplicates('ModelID')
    stranded = pd.Series(prof.reindex(rna.ProfileID).to_numpy(), index=rna.ModelID.to_numpy())
    res['by_stranded'] = {str(v): dict(n_lines=int(per['sample'].isin(stranded.index[stranded.eq(v)]).sum()),
                                       **primary(per[per['sample'].isin(stranded.index[stranded.eq(v)])], dev)['all'])
                          for v in sorted(stranded.dropna().unique())}
    res['family'] = 'Holm over 2: hd_dose vs naive, hd_dose vs infercnvpy (focal-loss recall, 5% FPR)'
    res['unit'] = 'cell line (DepMap ModelID), lines with >= 1 focal-loss window'
    res['score_version'] = SCORE_VERSION
    res['evidence_role'] = 'retrospective_regression; original frozen claims are not redefined'
    (out / 'primary.json').write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


def subarm_positions(inp):
    """Gene -> hg38 chromosome ('chrN') and start, from the groups' positions.tsv (one Ensembl 112 table)."""
    p = pd.concat(pd.read_csv(f, sep='\t') for f in inp.glob('*/positions.tsv')).drop_duplicates('gene')
    p['chromosome'] = 'chr' + p.chromosome.astype(str)
    return p.set_index('gene')[['chromosome', 'start']]


if __name__ == '__main__':
    main(*sys.argv[1:4])
