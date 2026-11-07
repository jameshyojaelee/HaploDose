#!/usr/bin/env python3
"""Merge per-specimen comparator runs (Numbat runs one specimen at a time) into one run layout for score.py.

From each run directory, the rows of its own specimen are kept from cells.tsv.gz, clones.tsv, arms.tsv and
differences.tsv (the other specimens' cells are NA there). Only completed runs (qc.json present) are merged; a
specimen without one is listed in OUT/missing_specimens.txt and contributes nothing (score.py then takes the patient
score over the specimens that ran). No value is changed.
Usage: merge_runs.py OUT SPECIMEN=RUN_DIR [SPECIMEN=RUN_DIR ...]
"""
import json
import sys
from pathlib import Path

import pandas as pd

TABLES = {'cells.tsv.gz': dict(compression='gzip'), 'clones.tsv': {}, 'arms.tsv': {}, 'differences.tsv': {}}


def main(out, pairs):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    merged = {t: [] for t in TABLES}
    used, missing = {}, []
    for p in pairs:
        spec, run = p.split('=', 1)
        run = Path(run)
        if not (run / 'qc.json').exists():
            missing.append(spec)
            continue
        used[spec] = str(run.resolve())
        for t in TABLES:
            f = run / t
            if f.exists() and f.stat().st_size and (t != 'differences.tsv' or f.read_text().strip()):
                d = pd.read_csv(f, sep='\t', dtype={'clone': str, 'clone_a': str, 'clone_b': str})
                merged[t].append(d[d.specimen.astype(str).eq(spec)])
    for t, kw in TABLES.items():
        if merged[t]:
            pd.concat(merged[t], ignore_index=True).to_csv(out / t, sep='\t', index=False, **kw)
    (out / 'missing_specimens.txt').write_text(''.join(f'{s}\n' for s in missing))
    (out / 'merged_from.json').write_text(json.dumps(used, indent=2) + '\n')
    if used:
        (out / 'qc.json').write_text(json.dumps(dict(merged_specimens=len(used), missing_specimens=missing)) + '\n')
    print(f'{len(used)} specimens merged, {len(missing)} missing')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])
