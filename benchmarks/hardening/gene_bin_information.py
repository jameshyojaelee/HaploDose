"""Matched gene/summed-gene/moment-NB oracles on original recovery failures."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd
from scipy import special, stats
from haplodose.provenance import start_run, finish_run
from haplodose.simulate import simulate


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def metrics(loglik, priors, truth):
    posterior = special.softmax(loglik+np.log(priors), axis=1)
    predicted = posterior.argmax(axis=1)
    normal, tumour = truth == 0, truth != 0
    score = special.log_softmax(loglik+np.log(priors), axis=1)
    return dict(all_correct=int(np.sum(predicted == truth)), all_n=len(truth),
        tumour_correct=int(np.sum(predicted[tumour] == truth[tumour])), tumour_n=int(tumour.sum()),
        normal_false_calls=int(np.sum(predicted[normal] != 0)), normal_n=int(normal.sum()),
        tumour_detected=int(np.sum(predicted[tumour] != 0)),
        mean_true_class_log_score=float(score[np.arange(len(truth)), truth].mean()))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--seed', type=int, choices=[20261108, 20261109, 20261110], required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    helper = load('positive_NB_sum', root/'gene_bin_reference.py')
    old = a.task.parent/'recovery_next_20261107'
    frozen = a.task.parent/f'recovery_20261107/recovery_panel_v1/seed_{a.seed}/data/counts.h5'
    mapping_path = old/f'trace_{a.seed}/fit/genes.tsv.gz'
    previous_cells = old/f'oracle_v1/cells_{a.seed}.tsv'
    rec = start_run(a.out, dict(script=Path(__file__), helper=root/'gene_bin_reference.py',
        generator=Path(simulate.__code__.co_filename), observations=frozen, mapping=mapping_path,
        previous_cells=previous_cells, previous_summary=old/'oracle_v1/summary.json',
        design=a.task/'DESIGN.md', numerical_gate=a.task/'numerics_v1/report.json'), dict(seed=a.seed),
        dict(scope='Known-generator measurement diagnostic only, biological n=0'), mode='gene_bin_information')
    try:
        assert json.loads((a.task/'numerics_v1/report.json').read_text())['all_checks_passed']
        captured = []

        def observer(frame, event, result):
            if event == 'return' and frame.f_code.co_name == 'draw' and frame.f_code.co_filename == simulate.__code__.co_filename:
                values = frame.f_locals
                if values.get('cols', {}).get('specimen') == 'S1':
                    captured.append(dict(label=values['cols']['truth_clone'], rate=values['rate'].copy(),
                        L=values['L'].copy(), gamma=values['gamma'].copy(), phi=values['phi'], depth=values['depth']))

        if sys.getprofile() is not None:
            raise RuntimeError('Cannot nest generator observation')
        sys.setprofile(observer)
        try:
            X, genes, meta, *_ = simulate(seed=a.seed)
        finally:
            sys.setprofile(None)
        with h5py.File(frozen) as source:
            for name, value in (('data', X.data), ('indices', X.indices), ('indptr', X.indptr)):
                assert np.array_equal(value, source['X/'+name][:]), name
            assert np.array_equal(meta.cell_id.to_numpy(), source['cell_ids'].asstr()[:])
            assert np.array_equal(np.asarray(genes), source['genes'].asstr()[:])
        assert [v['label'] for v in captured] == ['normal', 'clone1', 'clone2']
        rows = np.flatnonzero(meta.specimen.eq('S1'))
        counts = X[rows].toarray()
        mapping = pd.read_csv(mapping_path, sep='\t')
        mapping = mapping[mapping.bin >= 0]
        columns, bins = mapping.column.to_numpy(int), mapping.bin.to_numpy(int)
        assert len(set(columns)) == len(columns)
        assert np.array_equal(np.asarray(genes)[columns], mapping.gene.to_numpy())
        J = int(bins.max())+1
        assert np.array_equal(np.unique(bins), np.arange(J))
        L = np.concatenate([v['L'] for v in captured])
        phi, gamma, depth = captured[0]['phi'], captured[0]['gamma'], captured[0]['depth']
        assert all(v['phi'] == phi and v['depth'] == depth and np.array_equal(v['gamma'], gamma) for v in captured)
        priors = np.array([len(v['L']) for v in captured])/len(L)
        names = np.array(['normal', 'clone1', 'clone2'])
        truth = pd.Index(names).get_indexer(meta.iloc[rows].truth_clone)
        assert np.all(truth >= 0)
        observed = np.column_stack([counts[:, columns[bins == j]].sum(axis=1) for j in range(J)])
        gene_ll, exact_ll, moment_ll = (np.zeros((len(rows), 3)) for _ in range(3))
        gene_budget, exact_budget, moment_budget = (np.zeros_like(gene_ll) for _ in range(3))
        eps = np.finfo(float).eps
        for k, group in enumerate(captured):
            rates = group['rate'][None]*(L[:, None]/depth)**gamma[None]
            rates /= rates.sum(axis=1, keepdims=True)
            means = L[:, None]*rates
            pmf = stats.nbinom.logpmf(counts, 1/phi, 1/(1+phi*means))[:, columns]
            gene_ll[:, k] = pmf.sum(axis=1)
            operations = 128*len(columns)
            gene_budget[:, k] = operations*eps/(1-operations*eps)*np.maximum(1., np.abs(pmf).sum(axis=1))
            for j in range(J):
                selected = means[:, columns[bins == j]]
                value, budget = helper.summed_logpmf(observed[:, j], selected, np.full(selected.shape[1], 1/phi))
                exact_ll[:, k] += value
                exact_budget[:, k] += budget
                mean = selected.sum(axis=1)
                size = mean**2/(phi*(selected**2).sum(axis=1))
                value = stats.nbinom.logpmf(observed[:, j], size, size/(size+mean))
                moment_ll[:, k] += value
                moment_budget[:, k] += 128*selected.shape[1]*eps*np.maximum(1., np.abs(value))
            print(json.dumps(dict(seed=a.seed, class_label=group['label'], bins=J)), flush=True)
        loglik = dict(gene=gene_ll, summed=exact_ll, moment_NB=moment_ll)
        budgets = dict(gene=gene_budget, summed=exact_budget, moment_NB=moment_budget)
        summaries, predictions, margins = {}, {}, {}
        for method, value in loglik.items():
            score = value+np.log(priors)
            order = np.argsort(score, axis=1)
            margin = score[np.arange(len(rows)), order[:, -1]]-score[np.arange(len(rows)), order[:, -2]]
            numerical_budget = budgets[method].sum(axis=1)+128*J*eps*np.maximum(1., np.abs(value).sum(axis=1))
            if np.any(margin <= numerical_budget):
                raise ArithmeticError('Oracle decision overlaps numerical budget; independent recomputation required')
            summaries[method] = dict(**metrics(value, priors, truth), min_decision_margin=float(margin.min()),
                max_numerical_budget=float(numerical_budget.max()))
            predictions[method], margins[method] = score.argmax(axis=1), margin
        old_cells = pd.read_csv(previous_cells, sep='\t')
        assert np.array_equal(old_cells.cell_id.to_numpy(), meta.iloc[rows].cell_id.to_numpy())
        assert np.array_equal(names[predictions['gene']], old_cells.retained_genes_assignment.to_numpy())
        arrays = dict(cell_ids=meta.iloc[rows].cell_id.to_numpy(str), source_columns=columns, gene_bins=bins,
            retained_gene_ids=mapping.gene.to_numpy(str), source_rows=rows, priors=priors, truth=truth,
            L=L, gamma=gamma, phi=np.array(phi), depth=np.array(depth),
            class_rates=np.stack([v['rate'] for v in captured]), bin_counts=observed)
        for method in loglik:
            arrays[method+'_loglik'], arrays[method+'_budget'] = loglik[method], budgets[method]
            arrays[method+'_prediction'], arrays[method+'_margin'] = predictions[method], margins[method]
        np.savez_compressed(a.out/'evidence.npz', **arrays)
        pairs = []
        for first, second in (('gene', 'summed'), ('summed', 'moment_NB')):
            table = np.zeros((2, 2), int)
            np.add.at(table, ((predictions[first] == truth).astype(int), (predictions[second] == truth).astype(int)), 1)
            pairs.append(dict(first=first, second=second, correctness_table=table.tolist(),
                decision_disagreements=int(np.sum(predictions[first] != predictions[second]))))
        (a.out/'report.json').write_text(json.dumps(dict(seed=a.seed, exact_frozen_count_identity=True,
            retained_genes=len(columns), bins=J, summaries=summaries, paired=pairs,
            previous_retained_gene_decisions_identical=True, biological_specimens=0,
            production_changed=False, scientific_criteria_changed=False,
            scope='Known true generator parameters and classes; descriptive diagnostic, not full clone recovery or deployable performance'),
            indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU native frozen count replay, long-double positive PGF and SciPy NB')


if __name__ == '__main__':
    main()
