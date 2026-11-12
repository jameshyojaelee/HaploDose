"""Complete exposed production regression for the mandatory final-support repair."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose import clones, pipeline
from haplodose.io import write_json
from haplodose.provenance import file_record, start_run, finish_run


CONDITIONS = ('matched', 'state', 'two_clone')
REUSED = {('two_clone', 1), ('two_clone', 2), ('two_clone', 4), ('state', 3)}


def load(name):
    spec = importlib.util.spec_from_file_location('production_support_'+name,
        Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def location(task, condition, rep):
    seed = 202611098000+rep
    start = 0 if rep < 5 else 5
    original = task.parent/'pooled_prevalence_pilot_20261110'/(
        f'pilot_{condition}_{start:02d}_v1')/f'seed_{seed}'
    reuse = task/f'combined_phase8_{condition}_{rep}_v1/replay'/f'seed_{seed}'
    return original, reuse


def manifest_control(folder, original):
    manifest = json.loads((folder/'run_manifest.json').read_text())
    assert manifest['status'] == 'complete'
    for name in ('counts.h5', 'model_cells.tsv', 'positions.tsv'):
        expected = file_record(original/'input'/name)['sha256']
        matches = [r for r in manifest['input']['files'].values()
                   if Path(r['path']).name == name]
        assert len(matches) == 1 and matches[0]['sha256'] == expected, name
    return manifest


def feature_control(before, after):
    for name in ('genes.tsv.gz', 'bins.tsv'):
        pd.testing.assert_frame_equal(pd.read_csv(before/name, sep='\t'),
            pd.read_csv(after/name, sep='\t'), check_exact=True)
    np.testing.assert_array_equal(np.load(before/'bin_counts.npy'), np.load(after/'bin_counts.npy'))
    return dict(genes_bins_counts_exact=True)


def worker(args):
    assert args.condition in CONDITIONS and (args.start, args.stop) in ((0, 5), (5, 10))
    reader, observer = load('combined_saved_outcomes'), load('final_membership_full_recovery')
    assert reader.assignments([[324, 35, 1], [118, 100, 22]]) == 424
    inputs = dict(script=Path(__file__), design=args.task/'PRODUCTION_REGRESSION_DESIGN.md',
        scorer=Path(reader.__file__), observer=Path(observer.__file__),
        installed_wheel=args.task/'validation_23169639/installed_artifact.json',
        tests=args.task/'validation_23169639/pytest.txt')
    text = inputs['tests'].read_text()
    assert '498 passed' in text and 'failed' not in text.strip().splitlines()[-1]
    for rep in range(args.start, args.stop):
        original, reuse = location(args.task, args.condition, rep)
        for name in ('counts.h5', 'model_cells.tsv', 'truth_cells.tsv', 'positions.tsv'):
            inputs[f'{rep}_{name}'] = original/'input'/name
        inputs[f'{rep}_old_scores'] = original.parent/'results.json'
        for name in ('run_manifest.json', 'cells.tsv.gz', 'clones.tsv', 'differences.tsv',
                     'log.tsv', 'genes.tsv.gz', 'bins.tsv', 'bin_counts.npy'):
            inputs[f'{rep}_old_{name}'] = original/'production/fit'/name
        if (args.condition, rep) in REUSED:
            for path in (reuse/'production/fit').iterdir():
                if path.is_file():
                    inputs[f'{rep}_reuse_{path.name}'] = path
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, inputs, dict(condition=args.condition, start=args.start,
        stop=args.stop, minimum=30, phase=8, fresh_data=False),
        dict(scope='Complete exposed production regression; no new specimen-error guarantee'),
        mode='production_support_regression')
    try:
        rows = []
        for rep in range(args.start, args.stop):
            original, reuse = location(args.task, args.condition, rep)
            before = original/'production/fit'
            manifest_control(before, original)
            support = []
            reused = (args.condition, rep) in REUSED
            if reused:
                case = reuse
            else:
                case = args.out/f'seed_{202611098000+rep}'
                (case/'production').mkdir(parents=True)
                with observer.observe_support(support):
                    pipeline.run(original/'input/counts.h5', original/'input/model_cells.tsv',
                        original/'input/positions.tsv', case/'production/fit',
                        params=clones.CloneParams(split_mode='thin'), chunk=97,
                        counts_kind='raw', genome_build='hg38',
                        coordinate_system='one-based-closed', gene_namespace='symbol')
                assert support and all(r['all_supported'] and r['shapes_aligned'] for r in support)
            after = case/'production/fit'
            manifest_control(after, original)
            features = feature_control(before, after)
            truth = pd.read_csv(original/'input/truth_cells.tsv', sep='\t')
            old = reader.score(before, truth, args.condition)
            current = reader.score(after, truth, args.condition)
            official = next(r for r in json.loads((original.parent/'results.json').read_text())
                if r['replicate'] == rep and r['arm'] == 'production')
            assert all(old[k] == official[k] for k in old)
            preservation = reader.output_control(original, case)
            log = pd.read_csv(after/'log.tsv', sep='\t')
            removals = log.loc[log.action.eq('drop_small_final')].dropna(axis=1, how='all')
            if removals.empty:
                assert preservation['strict_original_preservation_passed'], (args.condition, rep)
            profiles = pd.read_csv(after/'clones.tsv', sep='\t')
            assert profiles.normal_like.dtype.kind == 'b'
            tumour = profiles.loc[~profiles.normal_like, 'n_eff']
            assert len(tumour) and np.isfinite(tumour).all() and (tumour >= 30).all()
            row = dict(condition=args.condition, replicate=rep, before=old, after=current,
                reused_complete_fit=reused, returned_support_checks=support,
                final_exported_support_passed=True, features=features,
                production_preservation=preservation, final_support_removals=removals.to_dict('records'),
                after_manifest=file_record(after/'run_manifest.json'))
            rows.append(row)
            write_json(rows, args.out/'cases.json')
            print(json.dumps(dict(condition=args.condition, replicate=rep, reused=reused,
                old_groups=old['reported_groups'], groups=current['reported_groups'],
                old_correct=old['assignment_correct'], correct=current['assignment_correct'])), flush=True)
        write_json(dict(cases=rows, count=5, fresh_data=False, scientific_criteria_changed=False,
            independent_evaluation=False, false_split_control_established=False,
            biological_specimens=0, promotion=False), args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU corrected production; fixed raw-input regression and exact reuse')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


def summary(args):
    inputs = dict(script=Path(__file__), design=args.task/'PRODUCTION_REGRESSION_DESIGN.md')
    for condition in CONDITIONS:
        for start in (0, 5):
            path = args.task/f'production_regression_{condition}_{start:02d}_v1'
            inputs[f'{condition}{start}_report'] = path/'report.json'
            inputs[f'{condition}{start}_manifest'] = path/'run_manifest.json'
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, inputs, dict(conditions=CONDITIONS, cases_per_condition=10),
        dict(scope='Complete paired regression counts; not new calibration or biological replication'),
        mode='production_support_regression_summary')
    try:
        rows = []
        for condition in CONDITIONS:
            for start in (0, 5):
                assert json.loads(inputs[f'{condition}{start}_manifest'].read_text())['status'] == 'complete'
                report = json.loads(inputs[f'{condition}{start}_report'].read_text())
                assert report['count'] == 5
                rows.extend(report['cases'])
        assert len(rows) == 30 and len({(r['condition'], r['replicate']) for r in rows}) == 30
        assert sum(r['reused_complete_fit'] for r in rows) == 4
        assert all(r['features']['genes_bins_counts_exact'] and r['final_exported_support_passed'] for r in rows)
        descriptions = []
        for condition in CONDITIONS:
            subset = [r for r in rows if r['condition'] == condition]
            assert sorted(r['replicate'] for r in subset) == list(range(10))
            for version in ('before', 'after'):
                descriptions.append(dict(condition=condition, version=version, cases=10,
                    exact_chr12_recovery=sum(r[version]['exact_chr12_recovery'] for r in subset),
                    additional_group_events=sum(r[version]['null_error'] for r in subset),
                    assignment_correct=sum(r[version]['assignment_correct'] for r in subset),
                    assignment_n=6000))
        write_json(dict(cases=rows, summaries=descriptions, complete_cases=30,
            new_inference_cases=26, reused_complete_fits=4, fresh_data=False,
            independent_evaluation=False, false_split_control_established=False,
            biological_specimens=0, promotion=False, scientific_criteria_changed=False,
            original_strict_gate_failures_preserved=True), args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU complete paired30-case regression summary')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--mode', choices=('worker', 'summary'), required=True)
    parser.add_argument('--condition', choices=CONDITIONS)
    parser.add_argument('--start', type=int)
    parser.add_argument('--stop', type=int)
    args = parser.parse_args()
    (worker if args.mode == 'worker' else summary)(args)


if __name__ == '__main__':
    main()
