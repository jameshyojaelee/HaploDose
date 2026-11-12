"""Observe production matched gene pools without changing any inference decision."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    (a.out/'arrays').mkdir()
    root = Path(__file__).parent
    runner = load('production_recovery', root/'gene_thinning_recovery.py')
    hook = load('production_gene_observer', root/'gene_veto_weights.py')
    previous = a.task.parent/'exact_profile_merge_20261109/production_v1/panel'
    rec = start_run(a.out, dict(script=Path(__file__), design=a.task/'DESIGN.md',
        previous=previous/'summary.json'), {},
        dict(scope='Five exposed production recovery cases; observation only, no oracle'), mode='production_gene_trace')
    calls, argv = [], sys.argv
    try:
        sys.argv = [str(runner.__file__), '--out', str(a.out/'replay'), '--task',
                    str(a.task.parent/'gene_thinning_20261109'), '--mode', 'control', '--stage', 'hard']
        with hook.observe_weights(calls, a.out/'arrays', apply=False, weight_units='effective'):
            runner.main()
        current = a.out/'replay/panel'
        old = json.loads((previous/'summary.json').read_text())
        new = json.loads((current/'summary.json').read_text())
        for before, after in zip(old['cases'], new['cases']):
            assert before['seed']==after['seed'] and before['criteria']==after['criteria']
            for name in ('cells.tsv.gz', 'clones.tsv', 'differences.tsv', 'segments.tsv', 'profiles.tsv.gz'):
                pd.testing.assert_frame_equal(pd.read_csv(current/f'seed_{after["seed"]}/fit'/name, sep='\t'),
                    pd.read_csv(previous/f'seed_{before["seed"]}/fit'/name, sep='\t'), check_exact=True)
            assert np.array_equal(np.load(current/f'seed_{after["seed"]}/fit/posteriors.npy'),
                                  np.load(previous/f'seed_{before["seed"]}/fit/posteriors.npy'), equal_nan=True)
        write_json(dict(all_public_tables_and_posteriors_identical=True, all30_criteria_identical=True,
                        calls=calls, production_changed=False, biological_specimens=0), a.out/'report.json')
        finish_run(a.out, rec, backend='CPU unchanged production; observational effective-unit comparison')
    except BaseException as error:
        write_json(calls, a.out/'partial_calls.json')
        finish_run(a.out, rec, error=error)
        raise
    finally:
        sys.argv = argv


if __name__ == '__main__':
    main()
