"""Frozen full regression of the two-case-qualified partition-average expert."""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from haplodose import gene_segments
from haplodose.io import write_json
from haplodose.provenance import file_record, start_run, finish_run


CHUNKS = ((0, 46), (46, 92), (92, 138), (138, 181))
NEW = ('partition_average_alone', 'partition_average_stack')


def load(name):
    spec = importlib.util.spec_from_file_location('partition_regression_'+name,
        Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def preflight(task, kernel):
    """Require actual successful gates and the same executed numerical source."""
    manifest = json.loads((task/'panel_v1/run_manifest.json').read_text())
    panel = json.loads((task/'panel_v1/report.json').read_text())
    numerics = json.loads((task/'numerics_v1/report.json').read_text())
    assert manifest['status'] == 'complete'
    assert panel['gate'] and panel['cases'] == 2 and panel['previous_scores_exact'] == 16
    assert numerics['independent_enumeration_passed'] and numerics['target_exclusion_exact']
    assert numerics['singleton'] and numerics['shift_and_scale_invariance']
    assert tuple(r['case'] for r in panel['paired']) == kernel.LOCKED
    assert all(r['selected'] and r['new_expert_weight'] > 0 for r in panel['RNA_decisions'])
    for row in panel['paired']:
        bound = panel['score_numerics'][row['case']]['conditional_projection_error_bound']
        assert math.isfinite(bound)
        assert all(row['partition_average_stack']-row[m] > bound
            for m in ('gene_laplace_stack', 'pseudobulk_sm5', 'pseudobulk_mdl'))
    for key, path in (('script', Path(kernel.__file__)), ('old_kernel', Path(gene_segments.__file__))):
        assert file_record(path)['sha256'] == manifest['input']['files'][key]['sha256']
    for record in manifest['output']['files'].values():
        assert file_record(record['path'])['sha256'] == record['sha256']


def common_inputs(task):
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    return dict(script=Path(__file__), kernel=Path(__file__).with_name('partition_averaged_stack.py'),
        design=task/'REGRESSION_DESIGN.md', panel=task/'panel_v1/report.json',
        panel_manifest=task/'panel_v1/run_manifest.json', numerics=task/'numerics_v1/report.json',
        identities=root/'dev/reference_weighted_pseudobulk_20261109/regression_v1/predictions_locked.json',
        old_scores=root/'dev/gene_laplace_profile_20261110/summary_v2/scores.tsv')


def unchanged(run):
    for record in run['input']['files'].values():
        assert file_record(record['path'])['sha256'] == record['sha256']


def worker(args):
    assert (args.start, args.stop) in CHUNKS
    kernel = load('partition_averaged_stack')
    preflight(args.task, kernel)
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    inputs = common_inputs(args.task)
    inputs['truth'] = root/'results/kinker_input_v2_20261101_1551/depmap_cn.tsv.gz'
    inputs['old_kernel'] = Path(gene_segments.__file__)
    identities = json.loads(inputs['identities'].read_text())['cases']
    assert len(identities) == 181 and len({r['case'] for r in identities}) == 181
    selected = identities[args.start:args.stop]
    assert len(selected) == args.stop-args.start
    for row in selected:
        i = row['prediction_index']
        chunk = max(start for start, _ in CHUNKS if start <= i)
        inputs[f'RNA{i}'] = root/f'dev/gene_laplace_profile_20261110/regression_{chunk:03d}_v2/predictions_case{i}.tsv.gz'
        inputs[f'bins{i}'] = root/f'dev/reference_weighted_pseudobulk_20261109/regression_v1/bins_case{i}.tsv'
    run = start_run(args.out, inputs, dict(start=args.start, stop=args.stop, new_experts=1),
        dict(scope='181 exposed regression lines; relative RNA; no independent evaluation or calibrated uncertainty'),
        mode='partition_average_regression')
    try:
        predictions, decisions = [], []
        for row in selected:
            i = row['prediction_index']
            genes = pd.read_csv(inputs[f'RNA{i}'], sep='\t', float_precision='round_trip')
            bins = pd.read_csv(inputs[f'bins{i}'], sep='\t')
            assert genes.gene.is_unique and genes.case.eq(row['case']).all()
            gb = genes.bin.to_numpy()
            candidate, standalone, decision, arrays = kernel.fit(
                genes.rna_log_ratio.to_numpy(), gb, bins.smoothing_block.to_numpy())
            baseline = arrays['baseline'] if arrays else candidate
            np.testing.assert_array_equal(baseline[gb], genes.gene_laplace_stack)
            if not arrays:
                assert not decision['selected'] and decision['reason'] == 'old_expert_support_unavailable'
                decision['candidate_profile_rounding'] = 0.
            assert np.isfinite(candidate).all() and np.isfinite(standalone).all()
            genes['partition_average_stack'] = candidate[gb]
            genes['partition_average_alone'] = standalone[gb]
            genes.to_csv(args.out/f'predictions_case{i}.tsv.gz', sep='\t', index=False)
            np.savez_compressed(args.out/f'fit_case{i}.npz', **arrays)
            write_json(decision, args.out/f'decision_case{i}.json')
            decisions.append(dict(case=row['case'], prediction_index=i,
                **{k: v for k, v in decision.items() if k not in ('fold_records', 'stacking')}))
            predictions.append(genes)
            print(json.dumps(dict(case=row['case'], completed=len(decisions), selected=decision['selected'])), flush=True)
        # Every RNA prediction in this worker is immutable before the first DNA read.
        write_json(dict(cases=decisions, DNA_parsed=False,
            outputs=[file_record(args.out/f"predictions_case{r['prediction_index']}.tsv.gz") for r in selected]),
            args.out/'RNA_LOCKED.json')
        cn = pd.read_csv(inputs['truth'], sep='\t', index_col=0)
        assert cn.index.is_unique and cn.columns.is_unique
        rows, score_numerics = [], {}
        for row, genes, decision in zip(selected, predictions, decisions):
            ref = pd.Series(row['reference_DNA_weights'])
            g = genes[genes.gene.isin(cn.index)].copy()
            g['truth'] = cn.loc[g.gene, row['query_DNA_line'][0]].to_numpy()-cn.loc[g.gene, ref.index].to_numpy()@ref.to_numpy()
            g['mb'] = g.chromosome+':'+((g.start.to_numpy()+g.end.to_numpy())//2//1000000).astype(str)
            mb = g.groupby('mb')[['truth']+list(kernel.OLD)+list(NEW)].mean()
            for method in kernel.OLD+NEW:
                finite = np.isfinite(mb[[method, 'truth']]).all().all() and mb[method].nunique() > 1
                rows.append(dict(case=row['case'], method=method,
                    pearson=float(pearsonr(mb[method], mb.truth).statistic) if finite else np.nan,
                    windows=len(mb), common_genes=len(g)))
            centred = mb.partition_average_stack.to_numpy()-mb.partition_average_stack.mean()
            norm = float(np.linalg.norm(centred))
            perturbation = math.sqrt(len(mb))*decision['candidate_profile_rounding']
            bound = 2*perturbation/(norm-perturbation) if norm > perturbation else np.inf
            score_numerics[row['case']] = dict(conditional_projection_error_bound=bound,
                centred_prediction_norm=norm, perturbation_norm_bound=perturbation,
                scope='Fixed fitted weights and exported-profile arithmetic only; not statistical uncertainty')
        scores = pd.DataFrame(rows)
        old = pd.read_csv(inputs['old_scores'], sep='\t', float_precision='round_trip')
        before = old[old.case.isin([r['case'] for r in selected])].reset_index(drop=True)
        pd.testing.assert_frame_equal(before, scores[scores.method.isin(kernel.OLD)].reset_index(drop=True), check_exact=True)
        scores.to_csv(args.out/'scores.tsv', sep='\t', index=False)
        matrix = scores.pivot(index='case', columns='method', values='pearson')
        assert np.isfinite(matrix['partition_average_stack']).all()
        write_json(dict(cases=len(selected), ordered_cases=[r['case'] for r in selected],
            previous_scores_exact=len(before), old_profiles_exact=len(selected),
            selected=sum(d['selected'] for d in decisions), means=matrix.mean().to_dict(),
            score_numerics=score_numerics, primary='partition_average_stack',
            standalone_is_ablation_only=True, independent_evaluation=False,
            biological_calibration=False, production_changed=False), args.out/'report.json')
        unchanged(run)
        finish_run(args.out, run, backend='CPU same exact interval costs and frozen all-partition averaging')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


def summary(args):
    kernel = load('partition_averaged_stack')
    preflight(args.task, kernel)
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    inputs = common_inputs(args.task)
    inputs.update(scorer=Path(__file__).with_name('statistics.py'),
        groups=root/'dev/rna_predictive_smoothing_20261109/characterization_v1/RNA_choices_locked.json')
    for start, _ in CHUNKS:
        path = args.task/f'regression_{start:03d}_v1'
        for name in ('scores.tsv', 'report.json', 'run_manifest.json', 'RNA_LOCKED.json'):
            inputs[f'{start}_{name}'] = path/name
    run = start_run(args.out, inputs, dict(family_size=8, seed=2026111003, n_boot=20000),
        dict(scope='Fixed eight conditional comparisons on exposed lines; eight reference groups, not independent validation'),
        mode='partition_average_summary')
    try:
        cases = [r['case'] for r in json.loads(inputs['identities'].read_text())['cases']]
        score_frames, counts, numerics = [], [], {}
        for start, stop in CHUNKS:
            manifest = json.loads(inputs[f'{start}_run_manifest.json'].read_text())
            report = json.loads(inputs[f'{start}_report.json'].read_text())
            assert manifest['status'] == 'complete'
            assert report['cases'] == stop-start and report['ordered_cases'] == cases[start:stop]
            assert report['old_profiles_exact'] == stop-start and report['previous_scores_exact'] == 8*(stop-start)
            for record in manifest['output']['files'].values():
                assert file_record(record['path'])['sha256'] == record['sha256']
            score_frames.append(pd.read_csv(inputs[f'{start}_scores.tsv'], sep='\t', float_precision='round_trip'))
            counts.append(report['selected'])
            numerics.update(report['score_numerics'])
        scores = pd.concat(score_frames, ignore_index=True)
        assert len(scores) == 181*10 and not scores.duplicated(['case', 'method']).any()
        assert set(scores.case) == set(cases) and set(scores.method) == set(kernel.OLD+NEW)
        old = pd.read_csv(inputs['old_scores'], sep='\t', float_precision='round_trip')
        pd.testing.assert_frame_equal(old, scores[scores.method.isin(kernel.OLD)].reset_index(drop=True), check_exact=True)
        matrix = scores.pivot(index='case', columns='method', values='pearson').loc[cases]
        assert np.isfinite(matrix.partition_average_stack).all()
        groups = {r['case']: r['dependence_group'] for r in json.loads(inputs['groups'].read_text())['widths']}
        labels = np.array([groups[c] for c in cases])
        assert len(np.unique(labels)) == 8
        stats = load('statistics')
        records = []
        for name in kernel.OLD:
            record = stats.paired_mean_comparison(matrix.partition_average_stack, matrix[name],
                clusters=labels, seed=2026111003, n_boot=20000, family_size=8)
            # Preserve an unavailable comparator in the family; do not shrink its population.
            if record['n_missing']:
                record.update(testable=False, p=1., reason='Incomplete fixed 181-line comparator population')
            record.update(a='partition_average_stack', b=name, endpoint='pearson', independent_evaluation=False)
            records.append(record)
        records = stats.adjust_family(records)
        assert set(numerics) == set(cases)
        arithmetic_resolved = all(math.isfinite(v['conditional_projection_error_bound']) for v in numerics.values())
        conditional_gate = bool(arithmetic_resolved and all(r['superior'] for r in records))
        scores.to_csv(args.out/'scores.tsv', sep='\t', index=False)
        matrix.to_csv(args.out/'paired_profiles.tsv', sep='\t')
        write_json(dict(cases=181, old_profiles_exact=181, previous_scores_exact=1448,
            means=matrix.mean().to_dict(), selected=sum(counts), family_size=8, contrasts=records,
            conditional_gate=conditional_gate, finite_export_arithmetic_bounds=arithmetic_resolved,
            score_numerics=numerics, primary='partition_average_stack', standalone_is_ablation_only=True,
            reference_groups=8, conditional_bootstrap=True, independent_evaluation=False,
            biological_calibration=False, default_promoted=False,
            decision='Permit separately verified opt-in package integration' if conditional_gate else
                'Reject promotion; preserve negative results, no parameter/comparator tuning'), args.out/'report.json')
        unchanged(run)
        finish_run(args.out, run, backend='CPU fixed eight paired means and whole-reference-group bootstrap')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--mode', choices=('worker', 'summary'), required=True)
    parser.add_argument('--start', type=int)
    parser.add_argument('--stop', type=int)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    if args.mode == 'worker':
        worker(args)
    else:
        summary(args)


if __name__ == '__main__':
    main()
