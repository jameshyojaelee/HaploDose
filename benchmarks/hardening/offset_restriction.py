"""Prespecified no-free-offset ablation; no query truth or package default change."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from haplodose.clones import SpecimenModel
from haplodose.provenance import start_run, finish_run


@contextmanager
def zero_offset(trace):
    original = SpecimenModel._update_offset

    def update(model, R, doses, keep):
        trace.append(dict(cells=model.rows.n, bins=model.J, components=len(doses),
            max_previous_offset=float(np.max(np.abs(model.offset)))))
        model.offset = np.zeros(model.J)

    SpecimenModel._update_offset = update
    try:
        yield
    finally:
        SpecimenModel._update_offset = original


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    path = Path(__file__).with_name('membership_stage_panel.py')
    spec = importlib.util.spec_from_file_location('recovery', path)
    recovery = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(recovery)
    rec = start_run(a.out, dict(script=Path(__file__), runner=path, design=a.task/'DESIGN.md'), {},
        dict(scope='Original five development cases; independent query normal labels unavailable'), mode='offset_restriction')
    trace, previous = [], sys.argv
    try:
        sys.argv = [str(path), '--out', str(a.out/'panel'), '--stage', 'hard',
            '--panel-root', str(a.task.parent/'recovery_20261107/recovery_panel_v1'),
            '--assertion-source', str(a.task.parent/'recovery_20261107/baseline/tests/test_split_mode.py')]
        with zero_offset(trace):
            recovery.main()
        (a.out/'offset_trace.json').write_text(json.dumps(trace, indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finally:
        sys.argv = previous
    finish_run(a.out, rec, backend='CPU original full workflow, shared offset restricted to zero')


if __name__ == '__main__':
    main()
