#!/usr/bin/env python3
"""How much the low-depth truth loses: SNU-601 low-depth clones (SRA cell-tagged subset, scdna_lowdepth.py)
vs full-depth clones (original BAM, scdna_cn.py). Cells matched by read-group ID (same IDs in both).

Reports: cells in both; clone counts; ARI on cells assigned to a clone in both; contingency table; arms that
differ between clones (clone_differences.tsv) in each, with overlap. The low-depth clone set is a coarser
partition, so arms it calls should be a subset of the full-depth arms; recall = shared / full arms.
Usage: compare_depth.py FULL_LINE_DIR LOW_LINE_DIR DEST
"""
import json
import sys
from pathlib import Path

import pandas as pd
from sklearn.metrics import adjusted_rand_score


def main(full, low, dest):
    full, low, dest = Path(full), Path(low), Path(dest)
    dest.mkdir(parents=True, exist_ok=False)
    f = pd.read_csv(full / 'cells.tsv', sep='\t').set_index('rg')
    l = pd.read_csv(low / 'cells.tsv', sep='\t').set_index('rg')
    both = f.index.intersection(l.index)
    j = pd.DataFrame({'full': f.loc[both, 'clone'], 'low': l.loc[both, 'clone']})
    a = j[~j.full.isin(['excluded', 'unassigned']) & ~j.low.isin(['excluded', 'unassigned'])]
    pd.crosstab(a.full, a.low).to_csv(dest / 'contingency.tsv', sep='\t')
    arms = lambda d: set(pd.read_csv(d / 'clone_differences.tsv', sep='\t').arm)
    fa, la = arms(full), arms(low)
    s = dict(cells_full=int(len(f)), cells_low=int(len(l)), cells_both=int(len(both)), cells_assigned_both=int(len(a)),
             clones_full=json.loads((full / 'qc.json').read_text())['clones'],
             clones_low=json.loads((low / 'qc.json').read_text())['clones'],
             ari=float(adjusted_rand_score(a.full, a.low)) if len(a) else None,
             diff_arms_full=sorted(fa), diff_arms_low=sorted(la), diff_arms_shared=sorted(fa & la),
             diff_arm_recall=len(fa & la) / len(fa) if fa else None,
             diff_arm_precision=len(fa & la) / len(la) if la else None)
    (dest / 'summary.json').write_text(json.dumps(s, indent=1) + '\n')
    print(json.dumps(s, indent=1))
    print(pd.crosstab(a.full, a.low).to_string())


if __name__ == '__main__':
    main(*sys.argv[1:])
