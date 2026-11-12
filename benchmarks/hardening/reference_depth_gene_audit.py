"""Scalar gene-counterfactual reconstruction and common-operator metric audit."""
import argparse
import hashlib
import importlib.util
import json
import math
import statistics
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.provenance import start_run, finish_run

spec = importlib.util.spec_from_file_location('scalar', Path(__file__).with_name('depth_reference_estimation_audit.py'))
scalar = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scalar)


def correlation(x, y):
    mx, my = statistics.mean(x), statistics.mean(y)
    vx = math.fsum((v-mx)**2 for v in x)
    vy = math.fsum((v-my)**2 for v in y)
    return math.fsum((a-mx)*(b-my) for a, b in zip(x, y))/math.sqrt(vx*vy)


def ranks(x):
    order = sorted(range(len(x)), key=lambda i: x[i])
    result = [0.]*len(x)
    start = 0
    while start < len(x):
        end = start+1
        while end < len(x) and x[order[end]] == x[order[start]]:
            end += 1
        for i in order[start:end]:
            result[i] = .5*(start+1+end)
        start = end
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    folder = a.task/'gene_fixture_v1'
    prior = a.task.parent/'profile_pool_export_20261109/parity_v1'
    original = a.task.parent/'profile_pool_mechanism_20261109/panel_v1'
    run = start_run(a.out, dict(script=Path(__file__), design=a.task/'GENE_FIXTURE_DESIGN.md',
        scalar=Path(scalar.__file__), predictions=folder/'count_predictions.npz',
        targets=folder/'truth_and_profiles.npz', metrics=folder/'metrics.tsv',
        source_means=original/'primary_means_and_profiles.npz'), {},
        dict(scope='One exposed RNA fixture, independent numerical audit'), mode='reference_depth_gene_audit')
    try:
        source = np.load(original/'primary_means_and_profiles.npz')
        M = source['gene_mean']
        Y = np.load(prior/'fitted_new/bin_counts.npy')
        meta = pd.read_csv(prior/'data/cells.tsv', sep='\t')
        cells = pd.read_csv(prior/'fitted_new/cells.tsv.gz', sep='\t')
        genes = pd.read_csv(prior/'fitted_new/genes.tsv.gz', sep='\t')
        ref = meta.role.eq('reference') & meta.ref_group.eq('type0')
        query = cells.specimen.eq('A') & cells.status.eq('assigned')
        saved = np.load(folder/'truth_and_profiles.npz')
        predictions = np.load(folder/'count_predictions.npz')
        selection = json.loads((folder/'selection.json').read_text())
        frozen = json.loads((folder/'PREDICTIONS_FROZEN.json').read_text())
        assert all(hashlib.sha256((folder/name).read_bytes()).hexdigest() == h for name, h in frozen.items())
        assert np.array_equal(ref, saved['reference_mask']) and np.array_equal(query, saved['query_mask'])
        assert np.array_equal(genes.bin, saved['gene_bin'])
        exposure = np.array([math.fsum(row) for row in M])
        idx = np.flatnonzero(ref)
        first = min(idx, key=lambda i: exposure[i])
        second = max(idx, key=lambda i: exposure[i])
        b = np.array([math.log(float(v)/float(M[first, -1])) for v in M[first]])
        g = np.array([(math.log(float(M[second, j])/float(M[second, -1]))-float(b[j]))
            /math.log(float(exposure[second])/float(exposure[first])) for j in range(M.shape[1])])
        origin = math.log(float(exposure[first]))
        budget = 512*M.size*np.finfo(float).eps
        assert np.max(np.abs(g-saved['true_gene_slope'])) <= budget
        assert np.max(np.abs(b-saved['true_gene_intercept'])) <= budget
        reference = scalar.scalar_rates(b, g, exposure[ref], origin)
        assert np.max(np.abs(reference/M[ref]-1)) <= budget
        counterfactual = scalar.scalar_rates(b, g, exposure[query], origin)
        assert np.max(np.abs(counterfactual/saved['counterfactual_gene']-1)) <= budget
        E = [math.fsum(float(counterfactual[i, j]) for i in range(len(counterfactual))
            for j in range(M.shape[1]) if int(genes.bin.iloc[j]) == bin_id) for bin_id in range(Y.shape[1])]
        Q = [math.fsum(float(M[i, j]) for i in np.flatnonzero(query)
            for j in range(M.shape[1]) if int(genes.bin.iloc[j]) == bin_id) for bin_id in range(Y.shape[1])]
        truth = [math.log(q/e) for q, e in zip(Q, E)]
        centre = statistics.median(truth)
        truth = [v-centre for v in truth]
        assert np.max(np.abs(np.array(truth)-saved['target'])) <= budget
        assert np.max(np.abs(np.array(E)/saved['counterfactual']-1)) <= budget
        regions = json.loads((original/'regions.json').read_text())
        metrics = pd.read_csv(folder/'metrics.tsv', sep='\t')
        assert len(metrics) == 8
        selected = min(selection['candidates'], key=lambda v: (v['score'], v['order']))['method']
        assert selected == selection['selected']
        for label in ('pooled', 'quintile', 'unpenalized', 'selected'):
            actual_label = selected if label == 'selected' else label
            Epred = predictions[actual_label+'__expected']
            if label != 'quintile':
                b, g = predictions[actual_label+'__b'], predictions[actual_label+'__gamma']
                independent = scalar.scalar_rates(b, g, Y[query].sum(axis=1), float(predictions['origin'])).sum(axis=0)
                assert np.max(np.abs(independent/Epred-1)) <= budget
            observed = [sum(int(Y[i, j]) for i in np.flatnonzero(query)) for j in range(Y.shape[1])]
            value = [math.log((y+.5)/(float(e)+.5)) for y, e in zip(observed, Epred)]
            centre = statistics.median(value)
            value = [v-centre for v in value]
            smooth = [None]*len(value)
            for start, end in regions:
                for j in range(start, end):
                    smooth[j] = math.fsum(value[start:end])/(end-start)
            assert None not in smooth
            for operator, prediction in [('none', value), ('shared_naive_boundaries', smooth)]:
                actual = metrics[(metrics.method == label) & (metrics.operator == operator)].iloc[0]
                mse = math.fsum((a-b)**2 for a, b in zip(prediction, truth))/len(truth)
                pearson = correlation(prediction, truth)
                spearman = correlation(ranks(prediction), ranks(truth))
                assert abs(actual.mse-mse) <= budget
                assert abs(actual.Pearson-pearson) <= budget
                assert abs(actual.Spearman-spearman) <= budget
        (a.out/'report.json').write_text(json.dumps(dict(gene_reference_means=60*M.shape[1],
            gene_query_counterfactual_means=160*M.shape[1], target_bins=Y.shape[1], scalar_metrics=24,
            biological_specimens=0, numerical_budget=budget, selected=selected,
            all_checks_passed=True, production_changed=False, scientific_criteria_changed=False), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, run, error=error)
        raise
    finish_run(a.out, run, backend='CPU scalar gene-law, pooling, ranks and metric references')


if __name__ == '__main__':
    main()
