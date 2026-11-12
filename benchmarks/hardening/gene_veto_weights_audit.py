"""Independent gene-pool sums and exact sign-test tails for unchanged decisions."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def exact_sign(ya, yb, gbin, mask, seg, minimum, alpha):
    informative = ya+yb >= 10
    if sum(informative) < minimum:
        return dict(passed=True, reason='insufficient_informative_genes')
    # Positive common library factors cancel when centering log ratios. Median
    # centering of an even number uses the geometric mean of the two middle ratios.
    ratio = [(Fraction(float(a))+Fraction(1,2))/(Fraction(float(b))+Fraction(1,2))
        for a,b in zip(ya,yb)]
    order = sorted(ratio[j] for j in np.flatnonzero(informative))
    lo, hi = order[(len(order)-1)//2], order[len(order)//2]
    selected = np.flatnonzero(informative & (gbin >= seg['first_bin']) & (gbin <= seg['last_bin'])
        & mask[np.clip(gbin, 0, len(mask)-1)])
    n = len(selected)
    k = sum(ratio[j]**2 > lo*hi if seg['direction'] == 'higher' else ratio[j]**2 < lo*hi for j in selected)
    tail = Fraction(sum(math.comb(n,j) for j in range(k,n+1)), 2**n)
    return dict(passed=bool(n < minimum or tail < Fraction(str(alpha))), n=n, k=k,
        tail=str(tail), reason='insufficient_segment_genes' if n < minimum else 'sign_test')


def sums(X, weights):
    # Scalar compensated sums, independent of the producer's sparse matrix product.
    return np.array([math.fsum(float(v)*float(weights[i]) for i,v in
        zip(X.indices[X.indptr[j]:X.indptr[j+1]], X.data[X.indptr[j]:X.indptr[j+1]]))
        for j in range(X.shape[1])])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--candidate', action='store_true')
    p.add_argument('--source', type=Path, help='Explicit immutable replay directory')
    p.add_argument('--production-changed', action='store_true')
    a = p.parse_args(); a.out.mkdir(exist_ok=False)
    source = a.source if a.source is not None else a.task/('candidate_v1' if a.candidate else 'diagnostic_v1')
    observed = json.loads((source/'report.json').read_text())
    inputs = dict(script=Path(__file__), design=a.task/'DESIGN.md', report=source/'report.json')
    for path in sorted(source.rglob('*')):
        if path.is_file() and path.suffix in ('.npz','.json','.tsv','.gz'):
            inputs[str(path.relative_to(source))] = path
    previous = a.task.parent/'gene_thinning_weighted_20261109/both_v1'
    inputs['previous_folds'] = previous/'trace.json'
    for row in json.loads(inputs['previous_folds'].read_text()):
        inputs[f'previous_fold:{row["model"]}'] = Path(row['path'])
    for seed in (20261031,20261107,20261108,20261109,20261110):
        inputs[f'{seed}:truth'] = a.task.parent/f'recovery_20261107/recovery_panel_v1/seed_{seed}/data/cells.tsv'
        for name in ('cells.tsv.gz','differences.tsv'):
            inputs[f'{seed}:previous:{name}'] = previous/f'panel/seed_{seed}/fit'/name
    rec = start_run(a.out, inputs, {}, dict(scope='Unchanged inference; hypothetical weighted-veto population diagnostic'),
        mode='gene_veto_weights_audit')
    try:
        assert observed.get('apply', False) == a.candidate
        assert observed['all_fit_tables_and30_criteria_identical'] == (not a.candidate)
        current_folds = json.loads((source/'replay/trace.json').read_text())
        previous_folds = json.loads(inputs['previous_folds'].read_text())
        assert len(current_folds) == len(previous_folds) == 15
        for current, old in zip(current_folds,previous_folds):
            for name in ('model','case_seed','cells','genes','bins','seed','share','oracle_gene_size','veto_calls'):
                assert current[name] == old[name]
            with np.load(current['path']) as x, np.load(old['path']) as y:
                assert set(x.files) == set(y.files)
                for name in x.files:
                    assert np.array_equal(x[name],y[name]), (current['model'],name)
        results, models = [], {}
        for call in observed['calls']:
            if call['model'] not in models:
                with np.load(call['model']) as z:
                    X = sparse.csr_matrix((z['data'],z['indices'],z['indptr']), shape=tuple(z['shape'])).tocsc()
                    models[call['model']] = (X,z['sub'].copy(),z['gbin'].copy(),z['mask'].copy(),int(z['seed']),z['rows'].copy())
            X, sub, gbin, mask, seed, rows = models[call['model']]
            with np.load(call['arrays']) as z:
                assert len(set(rows)) == len(rows) and len(set(sub)) == len(sub)
                labels = z['labels'][sub]
                early = min(sum(labels == call['a']),sum(labels == call['b'])) < 10
                assert bool(early) == call['early_cell_abstention']
                hard = [sums(X, labels == group) for group in (call['a'],call['b'])]
                assert np.array_equal(hard[0],z['hard_a']) and np.array_equal(hard[1],z['hard_b'])
                check = exact_sign(*hard, gbin, mask, call['seg'],call['min_genes'],call['alpha'])
                assert early or check['passed'] == call['original_passed']
                row = dict(seed=seed,index=call['index'],caller=call['caller'],seg=call['seg'],
                    original_passed=call['original_passed'],hard=check,weighted_available=call['weighted']['available'])
                if 'weights' in z:
                    W = z['weights']; assert W.shape == (len(rows),2)
                    weighted = [sums(X,W[sub,i]) for i in (0,1)]
                    for value,name in zip(weighted,('weighted_a','weighted_b')):
                        original = z[name]
                        budget = 8*X.shape[0]*np.finfo(float).eps*np.maximum(1.,abs(value))
                        assert np.all(abs(original-value) <= budget)
                    if call.get('weight_units','raw') == 'effective':
                        c=[sum(Fraction(float(w))**2 for w in W[sub,j])/sum(Fraction(float(w)) for w in W[sub,j]) for j in (0,1)]
                        unit=np.array([float(x) for x in c])
                        allowance=64*len(sub)*np.finfo(float).eps
                        assert np.all(abs(z['weight_scale']-unit)<=allowance*np.maximum(1.,abs(unit)))
                        weighted=[values/unit[j] for j,values in enumerate(weighted)]
                        for value,name in zip(weighted,('effective_a','effective_b')):
                            assert np.all(abs(value-z[name])<=allowance*np.maximum(1.,abs(value)))
                    else:
                        assert call.get('weight_units','raw') == 'raw' and 'weight_scale' not in z
                    wc = exact_sign(*weighted,gbin,mask,call['seg'],call['min_genes'],call['alpha'])
                    assert wc['passed'] == call['weighted']['passed']
                    if 'n' in wc:
                        assert wc['n'] == call['weighted']['n'] and wc['k'] == call['weighted']['k']
                    row.update(weighted=wc,changed=bool(wc['passed'] != call['original_passed']))
                if a.candidate:
                    chosen = row['weighted']['passed'] if not early and row['weighted_available'] else call['original_passed']
                    assert call['applied_passed'] == chosen
                results.append(row)
        summary = json.loads((source/'replay/panel/summary.json').read_text())
        for row in summary['cases']:
            seed = row['seed']; folder = source/f'replay/panel/seed_{seed}/fit'
            fit = pd.read_csv(folder/'cells.tsv.gz',sep='\t')
            truth = pd.read_csv(inputs[f'{seed}:truth'],sep='\t')
            joined = fit.merge(truth[['cell_id','truth_clone']],on='cell_id',validate='one_to_one')
            diff = pd.read_csv(folder/'differences.tsv',sep='\t')
            if not a.candidate:
                for name,value in (('cells.tsv.gz',fit),('differences.tsv',diff)):
                    pd.testing.assert_frame_equal(value,pd.read_csv(inputs[f'{seed}:previous:{name}'],sep='\t'),check_exact=True)
            s1 = joined[joined.specimen == 'S1']; normal = s1.truth_clone == 'normal'
            fpr = sum(s1.loc[normal,'malignant_call'])/sum(normal)
            tpr = sum(s1.loc[~normal,'malignant_call'])/sum(~normal)
            sub = s1[s1.truth_clone.isin(['clone1','clone2']) & s1.malignant_call]
            groups = sorted(set(sub.clone)); n = min(2,len(groups))
            table = [[sum((sub.truth_clone == pop)&(sub.clone == g)) for g in groups] for pop in ('clone1','clone2')]
            accuracy = max(sum(table[i][i] for i in range(n)),sum(table[1-i][i] for i in range(n)))/len(sub)
            s2 = joined[(joined.specimen == 'S2') & joined.malignant_call]
            d = diff[(diff.specimen == 'S1') & ~diff.clone_a.str.endswith('_c0') & ~diff.clone_b.str.endswith('_c0')]
            chrom = sorted(set(d.chromosome))
            expected = dict(normal_fpr=(fpr,fpr < .1),tumour_tpr=(tpr,tpr > .9),
                s2_one_group=(len(set(s2.clone)),len(set(s2.clone)) == 1),s1_two_groups=(len(groups),len(groups) == 2),
                assignment_accuracy=(accuracy,accuracy > .65),difference_chromosomes=(chrom,chrom == ['chr12']))
            for name,(value,passed) in expected.items():
                assert row['criteria'][name]['value'] == value and row['criteria'][name]['passed'] == bool(passed)
            assert row['original_assertion_passed'] == all(passed for _,passed in expected.values())
        write_json(dict(all_checks_passed=True, original_criteria_checked=30, exact_folds_checked=15, veto_calls=len(results),
            results=results, biological_specimens=0, production_changed=a.production_changed, scientific_criteria_changed=False,
            interpretation='Population sensitivity only; weighted sign-tail is not a calibrated new test'),a.out/'report.json')
    except BaseException as error:
        finish_run(a.out,rec,error=error);raise
    finish_run(a.out,rec,backend='CPU independent sparse scalar sums and exact binomial tails')


if __name__ == '__main__':
    main()
