"""Challenge frozen depth regression with gene-level generation and bin aggregation."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from haplodose.clones import reference_depth_slopes
from haplodose.provenance import start_run, finish_run


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


method = load('method', 'reference_depth_softmax.py')
anchor = load('anchor', 'depth_reference_anchor.py')


def centre(value):
    return value-np.median(value)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    root = a.task.parent
    prior = root/'profile_pool_export_20261109/parity_v1'
    fit = prior/'fitted_new'
    means_path = root/'profile_pool_mechanism_20261109/panel_v1/primary_means_and_profiles.npz'
    regions_path = root/'profile_pool_mechanism_20261109/panel_v1/regions.json'
    run = start_run(a.out, dict(script=Path(__file__), model=Path(method.__file__),
        design=a.task/'GENE_FIXTURE_DESIGN.md', cells=fit/'cells.tsv.gz', input_cells=prior/'data/cells.tsv',
        counts=fit/'bin_counts.npy', genes=fit/'genes.tsv.gz', means=means_path, regions=regions_path), {},
        dict(scope='Exposed gene-level RNA fixture, no DNA or biological validation'), mode='reference_depth_gene_fixture')
    try:
        meta = pd.read_csv(prior/'data/cells.tsv', sep='\t')
        cells = pd.read_csv(fit/'cells.tsv.gz', sep='\t')
        genes = pd.read_csv(fit/'genes.tsv.gz', sep='\t')
        Y = np.load(fit/'bin_counts.npy')
        assert cells.cell_id.tolist() == meta.cell_id.tolist()
        reference = meta.role.eq('reference') & meta.ref_group.eq('type0')
        fitmask = reference & meta.ref_split.eq('fit')
        holdmask = reference & meta.ref_split.eq('holdout')
        querymask = cells.specimen.eq('A') & cells.status.eq('assigned')
        assert (int(fitmask.sum()), int(holdmask.sum()), int(querymask.sum())) == (30, 30, 160)
        train, hold, query = Y[fitmask].astype(float), Y[holdmask], Y[querymask]
        L = train.sum(axis=1)
        Q = query.sum(axis=1).astype(float)
        origin = float(np.log(L).mean())
        x = np.log(L)-origin
        mu = train.sum(axis=0)/train.sum()
        curvature = float(np.sum(L*x*x)/L.sum()*(1-mu@mu)/len(mu))
        family = [('pooled', None)]+[(f'ridge_{k}', curvature*10.**k) for k in range(4, -5, -1)]+[('unpenalized', 0.)]
        parameters, scores = {}, []
        for order, (label, penalty) in enumerate(family):
            if penalty is None:
                b, g, record = np.log(mu), np.zeros(len(mu)), None
            else:
                b, g, _, record = method.fit(train, penalty, 1e-10)
            predicted = method.means(b, g, Q, origin).sum(axis=0)
            score = method.score(hold, method.means(b, g, hold.sum(axis=1).astype(float), origin))
            parameters[label+'__b'] = b
            parameters[label+'__gamma'] = g
            parameters[label+'__expected'] = predicted
            scores.append(dict(method=label, order=order, score=score, solver=record))
        chosen = min(scores, key=lambda v: (v['score'], v['order']))['method']
        gamma = reference_depth_slopes(train, L, [np.arange(len(train))], 30)
        b, record = anchor.calibrate(mu, gamma, L, origin, 1e-12)
        parameters['quintile__expected'] = method.means(b, gamma, Q, origin).sum(axis=0)
        parameters['origin'] = np.array(origin)
        np.savez_compressed(a.out/'count_predictions.npz', **parameters)
        (a.out/'selection.json').write_text(json.dumps(dict(selected=chosen, candidates=scores,
            quintile_solver=record, fit_ids=meta.loc[fitmask, 'cell_id'].tolist(),
            holdout_ids=meta.loc[holdmask, 'cell_id'].tolist(), query_ids=meta.loc[querymask, 'cell_id'].tolist()), indent=2)+'\n')
        frozen = {name: hashlib.sha256((a.out/name).read_bytes()).hexdigest()
            for name in ('count_predictions.npz', 'selection.json')}
        (a.out/'PREDICTIONS_FROZEN.json').write_text(json.dumps(frozen, indent=2)+'\n')
        # Generating means are accessed only for post-freeze truth construction.
        truth_source = np.load(means_path)
        M = truth_source['gene_mean']
        assert M.shape[0] == len(Y) and M.shape[1] == len(genes)
        assert np.array_equal(truth_source['gene_columns'], np.arange(len(genes)))
        assert (genes.bin.to_numpy() >= 0).all() and int(genes.bin.max())+1 == Y.shape[1]
        exposure = M.sum(axis=1)
        r = np.flatnonzero(reference)
        first, second = int(r[np.argmin(exposure[r])]), int(r[np.argmax(exposure[r])])
        logshare = np.log(M/exposure[:, None])
        dx = np.log(exposure[second])-np.log(exposure[first])
        slope = ((logshare[second]-logshare[second, -1])-(logshare[first]-logshare[first, -1]))/dx
        intercept = logshare[first]-logshare[first, -1]
        true_origin = float(np.log(exposure[first]))
        reconstructed = method.means(intercept, slope, exposure[r], true_origin)
        operations = 256*M.size
        budget = operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)
        reference_error = float(np.max(np.abs(np.log(reconstructed/M[r]))))
        assert reference_error <= budget
        counterfactual_gene = method.means(intercept, slope, exposure[querymask], true_origin)
        counterfactual = np.bincount(genes.bin, weights=counterfactual_gene.sum(axis=0), minlength=Y.shape[1])
        query_mean = np.bincount(genes.bin, weights=M[querymask].sum(axis=0), minlength=Y.shape[1])
        assert abs(counterfactual.sum()-exposure[querymask].sum()) <= budget*exposure[querymask].sum()
        target = centre(np.log(query_mean/counterfactual))
        observed = query.sum(axis=0)
        expected = dict(pooled=parameters['pooled__expected'], quintile=parameters['quintile__expected'],
            unpenalized=parameters['unpenalized__expected'], selected=parameters[chosen+'__expected'])
        regions = json.loads(regions_path.read_text())
        metrics, stored = [], dict(target=target, observed=observed, counterfactual=counterfactual,
            counterfactual_gene=counterfactual_gene, query_mean=query_mean, gene_mean=M,
            reference_mask=reference.to_numpy(), query_mask=querymask.to_numpy(), gene_bin=genes.bin.to_numpy(),
            true_gene_slope=slope, true_gene_intercept=intercept, true_origin=np.array(true_origin))
        for label, E in expected.items():
            value = centre(np.log((observed+.5)/(E+.5)))
            segmented = np.empty_like(value)
            for start, end in regions:
                segmented[start:end] = value[start:end].mean()
            for operator, v in [('none', value), ('shared_naive_boundaries', segmented)]:
                metrics.append(dict(method=label, operator=operator, mse=float(np.mean((v-target)**2)),
                    Pearson=float(np.corrcoef(v, target)[0, 1]), Spearman=float(spearmanr(v, target).statistic)))
                stored[label+'__'+operator] = v
        pd.DataFrame(metrics).to_csv(a.out/'metrics.tsv', sep='\t', index=False)
        np.savez_compressed(a.out/'truth_and_profiles.npz', **stored)
        (a.out/'report.json').write_text(json.dumps(dict(selected=chosen, bins=Y.shape[1], biological_specimens=0,
            query_cells=160, fit_cells=30, holdout_cells=30, reference_reconstruction_log_error=reference_error,
            numerical_budget=budget, production_changed=False, scientific_criteria_changed=False,
            truth='Depth-standardized pooled RNA counterfactual, not DNA copies', metrics=metrics), indent=2)+'\n')
        assert all(hashlib.sha256((a.out/name).read_bytes()).hexdigest() == digest for name, digest in frozen.items())
    except BaseException as error:
        finish_run(a.out, run, error=error)
        raise
    finish_run(a.out, run, backend='CPU SciPy count regression and generating-gene counterfactual')


if __name__ == '__main__':
    main()
