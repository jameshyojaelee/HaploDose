"""Matched width ablations and reference-cluster paired uncertainty, regression only."""
import argparse
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom, pearsonr

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    base = root/'dev/reference_weighted_pseudobulk_20261109/regression_v1'
    current = args.task/'regression_v1'
    inputs = dict(script=Path(__file__), design=args.task/'CHARACTERIZATION_DESIGN.md',
        selector=Path(__file__).with_name('rna_predictive_smoothing.py'),
        identities=base/'predictions_locked.json', selected=current/'RNA_decisions_locked.json',
        prior_scores=current/'scores.tsv', truth=root/'results/kinker_input_v2_20261101_1551/depmap_cn.tsv.gz')
    identities = json.loads(inputs['identities'].read_text())['cases']
    decisions = {r['case']: r for r in json.loads(inputs['selected'].read_text())['cases']}
    for r in identities:
        i = r['prediction_index']
        inputs[f'RNA{i}'] = current/f'predictions_case{i}.tsv.gz'
        inputs[f'bins{i}'] = base/f'bins_case{i}.tsv'
    rec = start_run(args.out, inputs, dict(bootstrap_draws=20000, seed=2026110911,
        bootstrap_unit='connected component of shared query/reference biological line identifiers'),
        dict(scope='Post-hoc characterization of frozen regression; percentile uncertainty conditional on fixed features/references and independent clusters'),
        mode='rna_smoothing_characterization')
    try:
        spec = importlib.util.spec_from_file_location('frozen_selector', inputs['selector'])
        selector = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(selector)
        cases = [r['case'] for r in identities]
        assert len(cases) == len(set(cases)) == 181
        # Dependence groups use biological identifiers only, without DNA measurements.
        parent = list(range(len(cases)))
        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        owner = {}
        for j, r in enumerate(identities):
            for line in set(r['reference_lines']+r['query_DNA_line']):
                if line in owner:
                    a, b = find(j), find(owner[line])
                    parent[max(a, b)] = min(a, b)
                else:
                    owner[line] = j
        groups = np.array([find(j) for j in range(len(cases))])
        widths = sorted(set.intersection(*(set(r['widths']) for r in decisions.values())))
        assert 5 in widths
        risks = np.array([[decisions[c]['mean_prediction_error'][decisions[c]['widths'].index(w)] for w in widths] for c in cases])
        def choose_width(risk):
            gamma = (4*risks.size)*np.finfo(float).eps
            guard = gamma/(1-gamma)*max(1., float(np.max(np.abs(risk))))
            return widths[int(np.flatnonzero(risk <= risk.min()+2*guard)[0])]
        global_width = choose_width(risks.mean(axis=0))
        methods = ['adaptive_RNA']+[f'fixed_{w}' for w in widths]+['global_RNA_width']
        if len(set(groups)) > 1:
            methods += ['leave_reference_group_out_RNA_width']
        predictions, width_records = {}, []
        for j, r in enumerate(identities):
            i = r['prediction_index']
            g = pd.read_csv(inputs[f'RNA{i}'], sep='\t', float_precision='round_trip')
            bins = pd.read_csv(inputs[f'bins{i}'], sep='\t')
            track = np.array([np.median(g.loc[g.bin.eq(b), 'rna_log_ratio']) for b in bins.bin])
            blocks = bins.smoothing_block.to_numpy()
            fixed = {w: selector.running_median(track, blocks, w) for w in widths}
            assert np.array_equal(fixed[5][g.bin.to_numpy()], g.pseudobulk_sm5)
            selected_width = decisions[r['case']]['selected_width']
            assert np.array_equal(fixed[selected_width][g.bin.to_numpy()], g.predictive_smoothed)
            prediction = [fixed[selected_width]]+[fixed[w] for w in widths]+[fixed[global_width]]
            row = dict(case=r['case'], dependence_group=int(groups[j]), global_RNA_width=global_width)
            if len(set(groups)) > 1:
                leave_out = choose_width(risks[groups != groups[j]].mean(axis=0))
                prediction += [fixed[leave_out]]
                row['leave_reference_group_out_RNA_width'] = leave_out
            predictions[str(i)] = np.array(prediction)
            width_records.append(row)
        np.savez_compressed(args.out/'RNA_bin_predictions.npz', **predictions)
        write_json(dict(methods=methods, widths=width_records,
            predictions=file_record(args.out/'RNA_bin_predictions.npz')), args.out/'RNA_choices_locked.json')
        cn = pd.read_csv(inputs['truth'], sep='\t', index_col=0)
        scores = []
        for j, r in enumerate(identities):
            i = r['prediction_index']
            gene = pd.read_csv(inputs[f'RNA{i}'], sep='\t', float_precision='round_trip')
            g = gene.loc[gene.gene.isin(cn.index)].copy()
            weights = pd.Series(r['reference_DNA_weights'])
            g['truth'] = cn.loc[g.gene, r['query_DNA_line'][0]].to_numpy()-cn.loc[g.gene, weights.index].to_numpy()@weights.to_numpy()
            for k, method in enumerate(methods):
                g[method] = predictions[str(i)][k, g.bin.to_numpy()]
            assert np.isfinite(g[['truth']+methods].to_numpy()).all()
            g['mb'] = g.chromosome+':'+((g.start.to_numpy()+g.end.to_numpy())//2//1_000_000).astype(str)
            mb = g.groupby('mb')[['truth']+methods].mean()
            for method in methods:
                scores.append(dict(case=r['case'], dependence_group=int(groups[j]), method=method,
                    pearson=float(pearsonr(mb[method], mb.truth).statistic)))
        scores = pd.DataFrame(scores)
        scores.to_csv(args.out/'scores.tsv', sep='\t', index=False)
        matrix = scores.pivot(index='case', columns='method', values='pearson').loc[cases]
        old = pd.read_csv(inputs['prior_scores'], sep='\t', float_precision='round_trip').pivot(index='case', columns='method', values='pearson').loc[cases]
        assert np.array_equal(matrix.fixed_5, old.pseudobulk_sm5)
        assert np.array_equal(matrix.adaptive_RNA, old.predictive_smoothed)
        effects = []
        for method in methods:
            if method == 'adaptive_RNA':
                continue
            delta = matrix.adaptive_RNA-matrix[method]
            effects.append(dict(comparator=method, wins=int((delta > 0).sum()), ties=int((delta == 0).sum()),
                losses=int((delta < 0).sum()), mean_paired_change=float(delta.mean()), median_paired_change=float(delta.median())))
        uncertainty = dict(status='not estimable from one connected reference group')
        clusters = np.unique(groups)
        if len(clusters) > 1:
            delta = (matrix.adaptive_RNA-matrix.fixed_5).to_numpy()
            total = np.array([delta[groups == c].sum() for c in clusters])
            size = np.array([sum(groups == c) for c in clusters])
            rng = np.random.default_rng(2026110911)
            sampled = rng.integers(len(clusters), size=(20000, len(clusters)))
            bootstrap = total[sampled].sum(axis=1)/size[sampled].sum(axis=1)
            np.save(args.out/'cluster_bootstrap_mean_differences.npy', bootstrap)
            values = np.sort(bootstrap)
            interval = np.quantile(values, [.025, .975]).tolist()
            mc_intervals = []
            for probability in (.025, .975):
                lo = max(0, int(binom.ppf(.0125, len(values), probability))-1)
                hi = min(len(values)-1, int(binom.ppf(.9875, len(values), probability)))
                mc_intervals.append([float(values[lo]), float(values[hi])])
            uncertainty = dict(status='conditional percentile cluster bootstrap', independent_clusters_assumed=len(clusters),
                mean_paired_change=float(delta.mean()), percentile_95_interval=interval,
                quantile_Monte_Carlo_order_statistic_brackets=mc_intervals,
                interpretation='Descriptive regression uncertainty; fixed feature universe/reference design, no prospective superiority claim')
        write_json(dict(cases=181, reference_dependence_groups=len(clusters),
            group_sizes={str(c): int(sum(groups == c)) for c in clusters},
            means=scores.groupby('method').pearson.mean().to_dict(), effects=effects,
            global_width_selected_by_RNA=global_width, primary_vs_fixed5_uncertainty=uncertainty,
            all362_primary_baseline_scores_reproduced=True, independent_evaluation=False,
            production_changed=False), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU matched smoothing ablations and conditional reference-cluster bootstrap')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
