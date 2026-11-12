"""Exact original production fixture before using an additional compute node."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose import clones, pipeline
from haplodose.io import write_json
from haplodose.provenance import file_record, start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--version', choices=('historical', 'corrected'), required=True)
    args = parser.parse_args()
    root = args.task.parent.parent
    original = root/'dev/pooled_prevalence_pilot_20261110/pilot_two_clone_00_v1/seed_202611098002'
    expected = (original/'production/fit' if args.version == 'historical' else
        root/'dev/final_membership_support_20261110/combined_phase8_two_clone_2_v1/replay/seed_202611098002/production/fit')
    assert json.loads((expected/'run_manifest.json').read_text())['status'] == 'complete'
    inputs = dict(script=Path(__file__), design=args.task/'PORTABILITY_DESIGN.md',
        counts=original/'input/counts.h5', metadata=original/'input/model_cells.tsv',
        positions=original/'input/positions.tsv', expected_manifest=expected/'run_manifest.json')
    eligible = lambda p: p.suffix == '.npy' or p.name.endswith(('.tsv', '.tsv.gz'))
    names = sorted(p.name for p in expected.iterdir() if eligible(p))
    assert len(names) == 16
    inputs.update({name:expected/name for name in names})
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, inputs, dict(version=args.version, equality='exact',
        fixture_seed=202611098002, prospective_outcomes_read=False),
        dict(scope='One fixed exposed fixture per version; node portability only'),
        mode='production_support_portability')
    try:
        actual = args.out/'fit'
        pipeline.run(inputs['counts'], inputs['metadata'], inputs['positions'], actual,
            params=clones.CloneParams(split_mode='thin'), chunk=97, counts_kind='raw',
            genome_build='hg38', coordinate_system='one-based-closed', gene_namespace='symbol')
        assert json.loads((actual/'run_manifest.json').read_text())['status'] == 'complete'
        assert names == sorted(p.name for p in actual.iterdir() if eligible(p))
        for name in names:
            if name.endswith('.npy'):
                np.testing.assert_array_equal(np.load(expected/name), np.load(actual/name))
            else:
                pd.testing.assert_frame_equal(pd.read_csv(expected/name, sep='\t'),
                    pd.read_csv(actual/name, sep='\t'), check_exact=True)
        write_json(dict(passed=True, arrays_exact=2, tables_including_logs_exact=14,
            version=args.version, prospective_outcomes_read=False,
            scope='This executed fixed fixture;not universal cross-node numerical parity'),
            args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU installed frozen production;exact portability fixture')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
