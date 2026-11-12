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
           'infercnvpy', 'copykat', 'scevan', 'stacked_smoothed')
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
    research = root/'dev/rna_stacked_smoothing_20261110'
    base = root/'dev/reference_weighted_pseudobulk_20261109/regression_v1'
    inp = root/'results/kinker_input_v2_20261101_1551'
    inputs = dict(script=Path(__file__), design=args.task/'FIELD_DESIGN.md',
        identities=base/'predictions_locked.json', truth=inp/'depmap_cn.tsv.gz',
        centromeres=Path(__file__).parents[2]/'python/haplodose/data/centromeres_hg38.tsv',
        scorer=Path(__file__).parents[1]/'kinker/evaluate.py',
        statistics=Path(__file__).with_name('statistics.py'),
        dependence=root/'dev/rna_predictive_smoothing_20261109/characterization_v1/RNA_choices_locked.json')
    cases = json.loads(inputs['identities'].read_text())['cases']
    paths = {}
    for start in (0, 46, 92, 138):
        d = research/f'regression_{start:03d}_v3'
        for r in json.loads((d/'RNA_DECISIONS_LOCKED.json').read_text())['cases']:
            paths[r['case']] = d/f"predictions_case{r['prediction_index']}.tsv.gz"
    runs = {}
    for name, dirname in ROOTS.items():
        for folder in sorted((root/'results'/dirname).glob('r*_pool*')):
            if not folder.is_dir():
                continue
            available = [p for p in folder.iterdir() if p.is_dir() and (p/'genes.tsv.gz').is_file()]
            assert len(available) == 1, (name, folder, available)
            rotation, pool = folder.name.split('_', 1)
            runs[name, int(rotation[1:]), pool] = available[0]
        assert sum(k[0] == name for k in runs) == 16
    for key, directory in runs.items():
        for filename in ('genes.tsv.gz', 'cells.tsv.gz', 'profiles.tsv.gz'):
            inputs['_'.join(map(str, key))+'_'+filename] = directory/filename
    for rotation in (0, 1):
        directory = root/f'results/kinker_field_20261102_0529/haplodose_r{rotation}'
        for filename in ('genes.tsv.gz', 'cells.tsv.gz', 'profiles.tsv.gz'):
            inputs[f'haplodose{rotation}_{filename}'] = directory/filename
    for case, path in paths.items():
        inputs['RNA_'+case] = path
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, inputs, dict(family_size=22, n_boot=20000, seed=2026111002),
        dict(scope='Exposed retrospective RNA profiles; fixed references, eight shared biological groups; no independent validation'),
        mode='stacked_field_comparison')
    try:
        scorer = load('stacked_field_scorer', inputs['scorer'])
        stats = load('stacked_field_statistics', inputs['statistics'])
        coverage, outputs = [], []
        for r in cases:
            gene = pd.read_csv(paths[r['case']], sep='\t', float_precision='round_trip').set_index('gene')
            gene['gene_pseudobulk'] = gene.rna_log_ratio
            hd = root/f"results/kinker_field_20261102_0529/haplodose_r{r['rotation']}"
            used, meta, profiles = read_run(hd)
            is_query = meta.specimen.eq(r['case']) & meta.design.eq('single_line') & meta.role.eq('query')
            query = meta[is_query]
            assert len(query) == r['query_cells']
            assert query.dataset.nunique() == 1
            pool = query.dataset.iloc[0]
            is_ref = meta.dataset.eq(pool) & meta.role.eq('reference') & meta.ref_split.eq('fit')
            # The frozen ordered_cells record contains QUERY PLUS FIT REFERENCE,
            # not queries alone. Reconstruct its exact original role selection.
            assert meta.loc[is_query | is_ref, 'cell_id'].tolist() == r['ordered_cells']
            query_ids, ref_ids = query.cell_id.tolist(), meta.loc[is_ref, 'cell_id'].tolist()
            assert len(ref_ids) == r['reference_cells']
            selected = [('haplodose_gene', hd, 'log_ratio_gene')]+[
                (name, runs[name, r['rotation'], pool], 'log_ratio') for name in ROOTS]
            for name, directory, quantity in selected:
                used, meta, profiles = read_run(directory)
                query = meta[meta.specimen.eq(r['case']) & meta.design.eq('single_line') & meta.role.eq('query')]
                assert query.cell_id.tolist() == query_ids
                ref = meta[meta.dataset.eq(pool) & meta.role.eq('reference') & meta.ref_split.eq('fit')]
                assert ref.cell_id.tolist() == ref_ids and set(ref.line) == set(r['reference_lines'])
                dna_weights = ref.groupby('line').umi.sum()
                assert (dna_weights/dna_weights.sum()).to_dict() == r['reference_DNA_weights']
                # Fixed all-method gene intersection; a missing population/bin
                # stays missing, never reweights the surviving clone population.
                common = gene.index.intersection(used.gene, sort=False)
                gene = gene.loc[common].copy()
                use = used.set_index('gene').loc[common]
                prof, fraction, status = scorer.specimen_profile(
                    profiles[profiles.specimen.eq(r['case'])], query.clone, use.bin.unique(), quantity)
                gene[name] = use.bin.map(prof).to_numpy()
                coverage.append(dict(case=r['case'], method=name, population=len(query),
                                     fraction=fraction, status=status, quantity=quantity))
            destination = args.out/f"RNA_case{r['prediction_index']}.tsv.gz"
            gene[list(('chromosome', 'start', 'end')+METHODS)].to_csv(destination, sep='\t')
            outputs.append(destination)
        pd.DataFrame(coverage).to_csv(args.out/'population_coverage.tsv', sep='\t', index=False)
        write_json(dict(methods=METHODS, decisions='All predictions fixed; common genes by identity only',
                        outputs=[file_record(p) for p in outputs]), args.out/'RNA_LOCKED.json')
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
                record = stats.paired_mean_comparison(matrix.stacked_smoothed[complete], matrix[method][complete],
                    clusters=labels[complete], seed=2026111002, n_boot=20000, family_size=22)
                record.update(a='stacked_smoothed', b=method, endpoint=endpoint,
                              all_family_complete_cases=int(complete.sum()), missing_family_cases=int((~complete).sum()))
                if not complete.all():
                    record['testable'] = False
                    record['p'] = np.nan
                    record['failure_reason'] = 'Incomplete fixed comparator family; no superiority verdict'
                contrasts.append(record)
        contrasts = stats.adjust_family(contrasts)
        write_json(dict(cases=len(cases), family_size=22, means=means, contrasts=contrasts,
            reference_groups=8, independent_evaluation=False, biological_calibration=False,
            frozen_inputs_unchanged=True, interpretation='Conditional retrospective comparison at common support; no universal superiority'),
            args.out/'report.json')
        finish_run(args.out, run, backend='CPU frozen field profiles and whole-reference-group bootstrap')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
