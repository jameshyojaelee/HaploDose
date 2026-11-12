"""Paired regression of the production population repair on frozen null/power inputs."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

from haplodose.provenance import start_run, finish_run


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--condition', choices=['matched_nb', 'correlated_state', 'two_clones'], required=True)
    p.add_argument('--start', type=int, choices=[0, 100], required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    old = a.task.parent/'dispersion_membership_20261109'
    runner = old/'tradeoff_source_v1/benchmarks/hardening/dispersion_tradeoff.py'
    spec = importlib.util.spec_from_file_location('frozen_tradeoff', runner)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rec = start_run(a.out, dict(script=Path(__file__), runner=runner, design=a.task/'TRADEOFF_DESIGN.md',
        recovery=a.task/'recovery_population_v1/summary.json'), {},
        dict(scope='Frozen-input paired regression; not independent evaluation or biological validation'),
        mode='offset_tradeoff')
    previous = sys.argv
    try:
        recovery = json.loads((a.task/'recovery_population_v1/summary.json').read_text())
        assert next(v for v in recovery['cases'] if v['seed'] == 20261031)['original_assertion_passed']
        sys.argv = [str(runner), '--out', str(a.out/'panel'), '--task', str(old), '--arm', 'MAP',
                    '--condition', a.condition, '--start', str(a.start)]
        module.main()
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finally:
        sys.argv = previous
    finish_run(a.out, rec, backend='CPU frozen complete binned workflow with installed population-repair wheel')


if __name__ == '__main__':
    main()
