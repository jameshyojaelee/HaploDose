"""Matched-population/support comparison of frozen RNA stacks and field tools."""
import argparse
from functools import lru_cache
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


METHODS = ('pseudobulk', 'pseudobulk_sm5', 'pseudobulk_mdl', 'predictive_smoothed',
           'predictive_partition', 'gene_pseudobulk', 'haplodose_gene', 'infercnv_r',
           'infercnvpy', 'copykat', 'scevan', 'stacked_smoothed', 'gene_laplace_alone', 'gene_laplace_stack')
ROOTS = dict(infercnv_r='kinker_infercnv_r_20261101_2110',
             infercnvpy='kinker_infercnvpy_20261101_1552_v2',
             copykat='kinker_copykat_20261101_2246', scevan='kinker_scevan_20261101_2246')


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@lru_cache(None)
def read_run(path):
    values = [pd.read_csv(path/name, sep='\t', float_precision='round_trip')
              for name in ('genes.tsv.gz', 'cells.tsv.gz', 'profiles.tsv.gz')]
    assert values[0].gene.is_unique and values[1].cell_id.is_unique
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    prior = root/'dev/native_stacked_profile_20261110/field_v2'
    base = root/'dev/reference_weighted_pseudobulk_20261109/regression_v1'
    inputs = dict(script=Path(__file__), design=args.task/'FIELD_DESIGN.md',
        identities=base/'predictions_locked.json',
        truth=root/'results/kinker_input_v2_20261101_1551/depmap_cn.tsv.gz',
        centromeres=Path(__file__).parents[2]/'python/haplodose/data/centromeres_hg38.tsv',
        scorer=Path(__file__).parents[1]/'kinker/evaluate.py',
        statistics=Path(__file__).with_name('statistics.py'), previous_scores=prior/'scores.tsv',
        population_coverage=prior/'population_coverage.tsv', previous_lock=prior/'RNA_LOCKED.json',
        dependence=root/'dev/rna_predictive_smoothing_20261109/characterization_v1/RNA_choices_locked.json')
    cases = json.loads(inputs['identities'].read_text())['cases']
    for r in cases:
        i = r['prediction_index']
        chunk = max(start for start in (0,46,92,138) if start <= i)
        inputs[f'field{i}'] = prior/f'RNA_case{i}.tsv.gz'
        inputs[f'new{i}'] = args.task/f'regression_{chunk:03d}_v2/predictions_case{i}.tsv.gz'
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, inputs, dict(family_size=26, n_boot=20000, seed=2026111004),
        dict(scope='Former evaluation regression-only; frozen common support and missing-family rules'),
        mode='gene_laplace_field')
    try:
        scorer = load('gene_field_scorer', inputs['scorer'])
        stats = load('gene_field_statistics', inputs['statistics'])
        outputs = []
        for r in cases:
            i = r['prediction_index']
            gene = pd.read_csv(inputs[f'field{i}'], sep='\t', index_col=0, float_precision='round_trip')
            new = pd.read_csv(inputs[f'new{i}'], sep='\t', float_precision='round_trip').set_index('gene')
            assert gene.index.is_unique and new.index.is_unique
            np.testing.assert_array_equal(gene.stacked_smoothed, new.loc[gene.index, 'stacked_smoothed'])
            for name in METHODS[-2:]:
                gene[name] = new.loc[gene.index, name]
            path = args.out/f'RNA_case{i}.tsv.gz'
            gene.to_csv(path, sep='\t')
            outputs.append(path)
        write_json(dict(methods=METHODS, outputs=[file_record(p) for p in outputs],
                        prior_population_coverage=file_record(inputs['population_coverage'])), args.out/'RNA_LOCKED.json')
        # No DNA was used to fit, select, change or normalize any prediction.
        cn = pd.read_csv(inputs['truth'], sep='\t', index_col=0, float_precision='round_trip')
        assert cn.index.is_unique and cn.columns.is_unique
        cen = pd.read_csv(inputs['centromeres'], sep='\t', comment='#').set_index('chromosome').p_arm_end
        rows = []
        for r, path in zip(cases, outputs):
            gene = pd.read_csv(path, sep='\t', index_col=0, float_precision='round_trip')
            gene = gene.loc[gene.index.intersection(cn.index, sort=False)].copy()
            weights = pd.Series(r['reference_DNA_weights'])
            gene['truth'] = (cn.loc[gene.index, r['query_DNA_line'][0]].to_numpy()-
                             cn.loc[gene.index, weights.index].to_numpy()@weights.to_numpy())
            gene['arm'] = gene.chromosome+np.where(gene.start.to_numpy() < gene.chromosome.map(cen).to_numpy(), 'p', 'q')
            assert gene.chromosome.map(cen).notna().all()
            gene['mb'] = gene.chromosome+':'+((gene.start.to_numpy()+gene.end.to_numpy())//2//1_000_000).astype(str)
            arm = gene.groupby('arm')[['truth']+list(METHODS)].median()
            arm = arm.loc[gene.groupby('arm').size() >= scorer.MIN_GENES]
            mb = gene.groupby('mb')[['truth']+list(METHODS)].mean()
            for method in METHODS:
                for endpoint, table, function in [('pearson_1mb', mb, pearsonr), ('spearman_arm', arm, spearmanr)]:
                    finite = (len(table) >= 2 and np.isfinite(gene[[method, 'truth']]).all().all()
                              and np.isfinite(table[[method, 'truth']]).all().all())
                    nonconstant = finite and table[method].nunique() > 1 and table.truth.nunique() > 1
                    score = float(function(table[method], table.truth).statistic) if nonconstant else np.nan
                    rows.append(dict(case=r['case'], method=method, endpoint=endpoint, score=score,
                                     genes=len(gene), regions=len(table), unavailable=not nonconstant))
        scores = pd.DataFrame(rows)
        scores.to_csv(args.out/'scores.tsv', sep='\t', index=False)
        previous = pd.read_csv(inputs['previous_scores'], sep='\t', float_precision='round_trip')
        actual = pd.read_csv(args.out/'scores.tsv', sep='\t', float_precision='round_trip')
        pd.testing.assert_frame_equal(previous, actual[actual.method.isin(METHODS[:-2])].reset_index(drop=True), check_exact=True)
        assert len(previous) == 4344
        groups = {r['case']:r['dependence_group'] for r in json.loads(inputs['dependence'].read_text())['widths']}
        labels = np.array([groups[r['case']] for r in cases])
        assert len(np.unique(labels)) == 8
        contrasts, means = [], {}
        for endpoint in ('pearson_1mb', 'spearman_arm'):
            matrix = scores[scores.endpoint.eq(endpoint)].pivot(index='case', columns='method', values='score')
            matrix = matrix.loc[[r['case'] for r in cases], list(METHODS)]
            complete = np.isfinite(matrix).all(axis=1).to_numpy()
            means[endpoint] = matrix.loc[complete].mean().to_dict()
            for method in METHODS[:-1]:
                record = stats.paired_mean_comparison(matrix.gene_laplace_stack[complete], matrix[method][complete],
                    clusters=labels[complete], seed=2026111004, n_boot=20000, family_size=26)
                record.update(a='gene_laplace_stack', b=method, endpoint=endpoint,
                              all_family_complete_cases=int(complete.sum()), missing_family_cases=int((~complete).sum()))
                if not complete.all():
                    record['testable'] = False
                    record['p'] = np.nan
                    record['failure_reason'] = 'Incomplete fixed comparator family; no superiority verdict'
                contrasts.append(record)
        contrasts = stats.adjust_family(contrasts)
        write_json(dict(cases=len(cases), family_size=26, means=means, contrasts=contrasts,
            previous_scores_exact=4344, reference_groups=8, independent_evaluation=False, biological_calibration=False,
            frozen_inputs_unchanged=True, interpretation='Conditional retrospective comparison at common support; no universal superiority'),
            args.out/'report.json')
        finish_run(args.out, run, backend='CPU frozen field profiles and whole-reference-group bootstrap')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
