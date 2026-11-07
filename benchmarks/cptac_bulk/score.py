"""CPTAC-3 bulk sub-arm benchmark (prespecified claim 6, docs/STATUS.md): runs benchmarks/tcga_bulk/subarm.py
unchanged except for its DNA inputs, driven as benchmarks/ccle_bulk/score.py, then the prespecified primary test.

New executions use version-2 retrospective scoring. The historical version-1
claim described below is preserved in its frozen result directory.

DNA inputs swapped in (subarm.py otherwise identical: 10 Mb gene-anchored windows, focal = segment shorter than
25% of its arm with |TCN - ploidy| >= 1, DNA-neutral = |TCN - ploidy| < 0.5, thresholds at the 5th/95th
percentile of DNA-neutral windows per cohort):
  gene positions  hg38 (Ensembl 112, the positions.tsv table); GDC AscatNGS segments are GRCh38, so no liftover
  arms           hg38 UCSC cytoBand (as in ccle_bulk/score.py)
  segments       AscatNGS Copy_Number (total CN); ploidy = length-weighted mean total CN over chr1-22
Tie-break: random (seeded jitter), so every method sits at the same realised FPR.

Primary (claim 6): focal-loss recall at 5% realised FPR per tumour, hd_dose vs naive and hd_dose vs infercnvpy,
paired Wilcoxon signed-rank over tumours with >= 1 focal-loss window, Holm over these 2. Claim holds if both Holm
p < 0.05 and the mean paired difference is > 0. Bootstrap 95% CI of the mean difference (10,000 resamples of
tumours, seed 20261103). Descriptive only: per cohort, and excluding tumours whose RNA file is pooled.
Usage: score.py OUT INPUT_ROOT ICNV_ROOT
"""
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd

H = Path(__file__).resolve().parents[2]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ccle = _load('ccle_score', H / 'benchmarks/ccle_bulk/score.py')   # hg38_arms, primary(), subarm_positions
subarm = ccle.subarm


def main(out, inp, icnv):
    out, inp, icnv = Path(out), Path(inp), Path(icnv)
    seg = pd.read_csv(inp / 'absolute_segments.tsv.gz', sep='\t')
    tab = pd.read_csv(inp / 'absolute_ploidy.tsv', sep='\t', index_col=0)
    pos = ccle.subarm_positions(inp)
    subarm.hg19_genes = lambda: pos
    subarm.hg19_arms = ccle.hg38_arms
    subarm.absolute = lambda: (tab, seg)
    subarm.GTF37 = 'hg38 positions (Ensembl 112), no liftover; GDC CPTAC-3 AscatNGS segments (GRCh38)'
    groups = sorted(d.name for d in inp.iterdir() if (d / 'counts.tsv').exists())
    sys.argv = ['subarm.py', str(out), '--tiebreak', 'random', '--input-scale', 'raw_counts'] + [f'{t}={inp / t}:{icnv / t}' for t in groups]
    subarm.main()
    per = pd.read_csv(out / 'per_patient.tsv.gz', sep='\t')
    pooled = set(pd.read_csv(inp / 'pooled_cases.tsv', sep='\t').case)
    res = {}
    r = ccle.primary(per, pooled)     # 'all' and the same excluding the given set
    res['all'] = r['all']
    res['excluding_pooled_rna'] = r['excluding_dev_lines']
    res['by_cohort'] = {t: dict(n_tumours=int(per['type'].eq(t).sum()), **ccle.primary(per[per['type'].eq(t)], set())['all'])
                        for t in groups}
    res['family'] = 'Holm over 2: hd_dose vs naive, hd_dose vs infercnvpy (focal-loss recall, 5% FPR)'
    res['unit'] = 'tumour (CPTAC-3 case), tumours with >= 1 focal-loss window'
    res['score_version'] = ccle.SCORE_VERSION
    res['descriptive'] = 'All version-2 outputs are retrospective regression; historical claims remain frozen.'
    (out / 'primary.json').write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main(*sys.argv[1:4])
