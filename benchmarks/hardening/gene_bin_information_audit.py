"""Independent positive convolutions, raw identity and oracle scoring audit."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, special
from haplodose.io import CountSource
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def nb_formula(y, means, size):
    y, means, size = np.asarray(y, float), np.asarray(means, float), np.asarray(size, float)
    return (special.gammaln(y+size)-special.gammaln(size)-special.gammaln(y+1)
        -size*np.log1p(means/size)+y*(np.log(means/size)-np.log1p(means/size)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    helper = load('positive_convolution', root/'gene_bin_reference.py')
    inputs = dict(script=Path(__file__), helper=root/'gene_bin_reference.py', design=a.task/'DESIGN.md',
        gate=a.task/'numerics_v1/report.json')
    for seed in (20261108, 20261109, 20261110):
        folder = a.task/f'seed_{seed}_v1'
        inputs[f'{seed}:evidence'], inputs[f'{seed}:report'] = folder/'evidence.npz', folder/'report.json'
        source = a.task.parent/f'recovery_20261107/recovery_panel_v1/seed_{seed}/data'
        inputs[f'{seed}:counts'], inputs[f'{seed}:cells'] = source/'counts.h5', source/'cells.tsv'
        inputs[f'{seed}:mapping'] = a.task.parent/f'recovery_next_20261107/trace_{seed}/fit/genes.tsv.gz'
    rec = start_run(a.out, inputs, {}, dict(scope='Independent oracle mechanism audit, biological n=0'),
        mode='gene_bin_information_audit')
    try:
        assert json.loads(inputs['gate'].read_text())['all_checks_passed']
        cases, independent = [], []
        eps = np.finfo(float).eps
        for seed in (20261108, 20261109, 20261110):
            saved = np.load(inputs[f'{seed}:evidence'])
            report = json.loads(inputs[f'{seed}:report'].read_text())
            meta = pd.read_csv(inputs[f'{seed}:cells'], sep='\t')
            rows = np.flatnonzero(meta.specimen.eq('S1'))
            assert np.array_equal(rows, saved['source_rows'])
            assert np.array_equal(meta.iloc[rows].cell_id, saved['cell_ids'])
            truth = pd.Index(['normal', 'clone1', 'clone2']).get_indexer(meta.iloc[rows].truth_clone)
            assert np.array_equal(truth, saved['truth'])
            mapping = pd.read_csv(inputs[f'{seed}:mapping'], sep='\t')
            mapping = mapping[mapping.bin >= 0]
            assert np.array_equal(mapping.column, saved['source_columns'])
            assert np.array_equal(mapping.bin, saved['gene_bins'])
            with CountSource(inputs[f'{seed}:counts']) as source:
                counts = sparse.vstack([chunk for _, _, chunk in source.chunks(97)]).tocsr()[rows].toarray()
                assert np.array_equal(source.genes[saved['source_columns']], saved['retained_gene_ids'])
            columns, bins = saved['source_columns'], saved['gene_bins']
            J = int(bins.max())+1
            observed = np.column_stack([counts[:, columns[bins == j]].sum(axis=1) for j in range(J)])
            assert np.array_equal(observed, saved['bin_counts'])
            priors = np.bincount(truth, minlength=3)/len(truth)
            assert np.array_equal(priors, saved['priors'])
            selected = set()
            for label in range(3):
                members = np.flatnonzero(truth == label)
                selected.update([int(members[0]), int(members[-1])])
            for method in ('gene', 'summed', 'moment_NB'):
                selected.add(int(np.argmin(saved[method+'_margin'])))
            selected = np.array(sorted(selected))
            exact = np.zeros((len(selected), 3))
            reference_budgets = np.zeros_like(exact)
            for k in range(3):
                rates = saved['class_rates'][k][None]*(saved['L'][:, None]/float(saved['depth']))**saved['gamma'][None]
                means = saved['L'][:, None]*rates/rates.sum(axis=1, keepdims=True)
                phi = float(saved['phi'])
                value = nb_formula(counts, means, 1/phi)[:, columns].sum(axis=1)
                assert np.all(np.abs(value-saved['gene_loglik'][:, k]) <= saved['gene_budget'][:, k])
                moment = np.zeros(len(truth))
                for j in range(J):
                    group = means[:, columns[bins == j]]
                    mean = group.sum(axis=1)
                    variance = (group+phi*group**2).sum(axis=1)
                    size = mean**2/(variance-mean)
                    moment += nb_formula(observed[:, j], mean, size)
                    for index, cell in enumerate(selected):
                        coefficient = helper.convolved_logpmf(observed[cell, j], group[cell], np.full(group.shape[1], 1/phi))
                        exact[index, k] += coefficient
                        operations = 128*group.shape[1]*(int(observed[cell, j])+1)**2
                        reference_budgets[index, k] += operations*eps/(1-operations*eps)*max(1., abs(coefficient))
                moment_budget = 128*J*eps*np.maximum(1., np.abs(moment))
                assert np.all(np.abs(moment-saved['moment_NB_loglik'][:, k]) <= saved['moment_NB_budget'][:, k]+moment_budget)
            discrepancy = np.abs(exact-saved['summed_loglik'][selected])
            budget = reference_budgets+saved['summed_budget'][selected]
            assert np.all(discrepancy <= budget)
            assert np.array_equal((exact+np.log(priors)).argmax(axis=1), saved['summed_prediction'][selected])
            independent.append(dict(seed=seed, cell_ids=saved['cell_ids'][selected].tolist(),
                cells=len(selected), full_bin_class_convolutions=len(selected)*J*3,
                max_summed_loglik_error=float(discrepancy.max()), max_budget=float(budget.max())))
            predictions = {}
            for method in ('gene', 'summed', 'moment_NB'):
                ll = saved[method+'_loglik']
                score = ll+np.array([math.log(p) for p in priors])
                predicted = score.argmax(axis=1)
                assert np.array_equal(predicted, saved[method+'_prediction'])
                predictions[method] = predicted
                normal, tumour = truth == 0, truth != 0
                expected = dict(all_correct=int(np.sum(predicted == truth)), all_n=len(truth),
                    tumour_correct=int(np.sum(predicted[tumour] == truth[tumour])), tumour_n=int(tumour.sum()),
                    normal_false_calls=int(np.sum(predicted[normal] != 0)), normal_n=int(normal.sum()),
                    tumour_detected=int(np.sum(predicted[tumour] != 0)))
                for name, value in expected.items():
                    assert report['summaries'][method][name] == value
                logscore = []
                for row, label in zip(score, truth):
                    maximum = float(row.max())
                    logscore.append(float(row[label])-maximum-math.log(math.fsum(math.exp(float(v)-maximum) for v in row)))
                average = math.fsum(logscore)/len(logscore)
                assert abs(average-report['summaries'][method]['mean_true_class_log_score']) <= 128*len(truth)*eps*max(1., abs(average))
                margin = np.sort(score, axis=1)[:, -1]-np.sort(score, axis=1)[:, -2]
                assert np.array_equal(margin, saved[method+'_margin'])
                numerical_budget = saved[method+'_budget'].sum(axis=1)+128*J*eps*np.maximum(1., np.abs(ll).sum(axis=1))
                assert np.all(margin > numerical_budget)
            for row in report['paired']:
                table = [[int(np.sum(((predictions[row['first']] == truth) == bool(i)) &
                    ((predictions[row['second']] == truth) == bool(j)))) for j in range(2)] for i in range(2)]
                assert table == row['correctness_table']
                assert row['decision_disagreements'] == int(np.sum(predictions[row['first']] != predictions[row['second']]))
            cases.append(report)
            print(json.dumps(dict(seed=seed, audited=True, selected_cells=len(selected))), flush=True)
        (a.out/'report.json').write_text(json.dumps(dict(all_checks_passed=True, cases=cases,
            independent_convolutions=independent, biological_specimens=0, production_changed=False,
            scientific_criteria_changed=False,
            scope='Known-parameter raw measurement comparison, not independent biological/full-workflow validation'), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU independent raw counts, direct PMF convolutions and scalar scoring')


if __name__ == '__main__':
    main()
