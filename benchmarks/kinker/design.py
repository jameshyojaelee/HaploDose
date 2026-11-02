"""Rebuild the Kinker rotation design (prepare.write_design, with prepare.eligible) from an existing input
directory into a new one; counts.h5, positions.tsv, depmap_cn.tsv.gz and cells_assigned.tsv.gz are linked,
not copied. Usage: design.py INPUT_DIR NEW_DIR
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from prepare import eligible, write_design  # noqa: E402

LINKED = ['counts.h5', 'positions.tsv', 'depmap_cn.tsv.gz', 'cells_assigned.tsv.gz', 'assignment_check.tsv']


def main(inp, out):
    inp, out = Path(inp).resolve(), Path(out)
    out.mkdir(parents=True, exist_ok=False)
    for f in LINKED:
        (out / f).symlink_to(inp / f)
    check = pd.read_csv(inp / 'assignment_check.tsv', sep='\t', dtype={'pool': str})
    cells = pd.read_csv(inp / 'cells_assigned.tsv.gz', sep='\t', dtype={'pool': str})
    ok = eligible(check)
    check.assign(eligible=ok).to_csv(out / 'eligible_lines.tsv', sep='\t', index=False)
    print(f'eligible lines {int(ok.sum())} of {len(check)}; excluded:',
          check.loc[~ok, ['pool', 'line']].apply(lambda r: f'{r.pool}:{r.line}', axis=1).tolist())
    write_design(cells, check, out)


if __name__ == '__main__':
    main(*sys.argv[1:])
