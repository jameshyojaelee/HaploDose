"""Exact installed-wheel replay of the locked native-region bulk candidate."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.bulk import run_bulk
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    root = args.task.parents[1]
    inp = root/'dev/hardening_20261106_2119/profile_v6/raw_refit_000'
    native = root/'dev/native_predictive_profile_20261109/raw_bulk_v1/haplodose'
    candidate = root/'dev/programme_gene_profile_20261109/bulk_v1'
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md')
    for name in ('counts.tsv', 'samples.tsv', 'positions.tsv'):
        inputs['input:'+name] = inp/name
    for path in native.iterdir():
        if path.suffix == '.npy' or path.name.endswith(('.tsv', '.tsv.gz')):
            inputs['native:'+path.name] = path
    inputs['candidate_predictions'] = candidate/'RNA_predictions.npz'
    for i in range(68):
        inputs[f'candidate_context_{i:03d}'] = candidate/f'context_{i:03d}.npz'
    run = start_run(args.out, inputs, dict(default_change=False),
        dict(scope='Exact artifact replay; original cohort is regression-only'), mode='native_region_replay')
    try:
        out = run_bulk(inp/'counts.tsv', inp/'samples.tsv', inp/'positions.tsv', args.out/'haplodose',
            min_cpm=0., device='cpu', counts_kind='raw', genome_build='hg38',
            coordinate_system='one-based-closed', gene_namespace='symbol', rna_profile_smoothing='native-regions')
        controls = []
        for name, path in inputs.items():
            if not name.startswith('native:') or path.name in ('rna_log_ratio_predictive.npy',):
                continue
            if path.suffix == '.npy':
                np.testing.assert_array_equal(np.load(path), np.load(out/path.name))
            else:
                pd.testing.assert_frame_equal(pd.read_csv(path, sep='\t'),
                                              pd.read_csv(out/path.name, sep='\t'), check_exact=True)
            controls.append(path.name)
        with np.load(candidate/'RNA_predictions.npz') as z:
            np.testing.assert_array_equal(np.load(out/'rna_log_ratio_native_regions.npy'),
                                          z['programme_gene_native_regions'])
        genes = pd.read_csv(out/'genes.tsv.gz', sep='\t')
        samples = pd.read_csv(out/'samples.tsv', sep='\t', dtype={'sample':str})
        adjusted = np.load(out/'gene_rna_log_ratio_programme.npy')
        expected = np.load(out/'gene_fitted_expected.npy')
        assert adjusted.shape == expected.shape == (68, len(genes))
        for i in range(len(samples)):
            with np.load(candidate/f'context_{i:03d}.npz') as z:
                np.testing.assert_array_equal(adjusted[i], z['corrected_gene'])
                np.testing.assert_array_equal(expected[i], z['expected_gene'])
                np.testing.assert_array_equal(np.load(out/'rna_state_index.npy')[i], z['path'])
        metadata = json.loads((out/'run.json').read_text())
        assert metadata['rna_profile_smoothing']['calibration_status'] == 'unvalidated'
        write_json(dict(exact_native_controls=controls, samples=len(samples), genes=len(genes),
            candidate_region_profiles_exact=True, candidate_gene_effects_exact=True,
            candidate_gene_expectations_exact=True, independent_evaluation=False,
            default_changed=False, clone_decisions_changed=False, scientific_criteria_changed=False,
            outputs={name:file_record(out/name) for name in ('rna_log_ratio_native_regions.npy',
                'gene_rna_log_ratio_programme.npy', 'gene_fitted_expected.npy', 'rna_native_regions.json')}),
            args.out/'report.json')
        finish_run(args.out, run, backend='CPU installed noneditable wheel; exact native-region replay')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
