"""Joint-law counterexamples and known-parameter contexts for NB count thinning."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import special, stats
from haplodose.provenance import start_run, finish_run


def distribution(means, sizes, cutoff):
    support, result = np.arange(cutoff+1), np.array([1.])
    for mean, size in zip(means, sizes):
        component = stats.nbinom.pmf(support, size, size/(size+mean))
        result = np.convolve(result, component)[:cutoff+1]
    return result


def grid(means, sizes, split_size, cutoff):
    f = distribution(means, sizes, cutoff)
    index = np.arange(cutoff+1)
    total = index[:, None]+index[None]
    valid = total <= cutoff
    joint = np.where(valid, f[np.minimum(total, cutoff)]*
        stats.betabinom.pmf(index[:, None], total, split_size/2, split_size/2), 0.)
    marginal = distribution(np.asarray(means)/2, np.asarray(sizes)/2, cutoff)
    gene_joint = np.where(valid, marginal[:, None]*marginal[None], 0.)
    tail = float(sum(stats.nbinom.sf(cutoff//len(means), r, r/(r+m)) for m, r in zip(means, sizes)))
    numerical_budget = (64*(cutoff+1)**2+128*len(means)*(cutoff+1))*np.finfo(float).eps
    result = {}
    for label, matrix in (('bin', joint), ('gene', gene_joint)):
        product = matrix.sum(axis=1)[:, None]*matrix.sum(axis=0)[None]
        value = float(np.abs(matrix-product).sum()/2)
        error = 1.5*tail+numerical_budget
        result[label] = dict(TV_truncated=value, TV_lower=max(0., value-error),
            TV_upper=min(1., value+error), captured_mass=float(matrix.sum()))
    a, b = split_size/2, split_size/2
    quadrature = []
    for order in (64, 128):
        x, weights = special.roots_jacobi(order, b-1, a-1)
        share = (x+1)/2
        zero = np.prod([(r/(r+m*share))**r for m, r in zip(means, sizes)], axis=0)
        quadrature.append(float((weights/weights.sum())@zero))
    assert abs(quadrature[0]-quadrature[1]) <= numerical_budget
    assert abs(joint[0].sum()-quadrature[-1]) <= tail+numerical_budget
    return dict(cutoff=cutoff, omitted_mass_union_bound=tail, numerical_budget=numerical_budget,
        bin=result['bin'], gene=result['gene'], independent_zero_half_quadrature=quadrature,
        joint_low=joint[:2, :2].tolist(), gene_joint_low=gene_joint[:2, :2].tolist())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    prior = a.task.parent/'gene_bin_information_20261109'
    inputs = dict(script=Path(__file__), design=a.task/'DESIGN.md', previous_audit=prior/'audit_v1/report.json')
    for seed in (20261108, 20261109, 20261110):
        inputs[f'{seed}:parameters'] = prior/f'seed_{seed}_v1/evidence.npz'
        inputs[f'{seed}:mapping'] = a.task.parent/f'recovery_next_20261107/trace_{seed}/fit/genes.tsv.gz'
    rec = start_run(a.out, inputs, {}, dict(scope='Conditional gene NB joint law, not workflow error control'),
        mode='gene_thinning_law')
    try:
        assert json.loads(inputs['previous_audit'].read_text())['all_checks_passed']
        M, V, share = Fraction(5), Fraction(22), Fraction(1, 2)
        R = M*M/(V-M)
        joint = [[Fraction(1, 10), Fraction(13, 200)],
                 [Fraction(13, 200), Fraction(129, 1000)*2*share*(1-share)*R/(R+1)]]
        determinant = joint[0][0]*joint[1][1]-joint[0][1]*joint[1][0]
        covariance = share*(1-share)*R/(R+1)*(V+M*M-M)-share*(1-share)*M*M
        assert covariance == 0 and determinant != 0
        examples = []
        for label, means, sizes in (('heterogeneous', [1., 4.], [1., 1.]), ('NB_closure', [1., 1.], [1., 1.])):
            means, sizes = np.array(means), np.array(sizes)
            split_size = float(means.sum()**2/np.sum(means**2/sizes))
            cutoff = 2
            while sum(stats.nbinom.sf(cutoff//len(means), r, r/(r+m)) for m, r in zip(means, sizes)) > np.sqrt(np.finfo(float).eps):
                cutoff *= 2
            coarse, fine = grid(means, sizes, split_size, cutoff), grid(means, sizes, split_size, 2*cutoff)
            assert abs(coarse['bin']['TV_truncated']-fine['bin']['TV_truncated']) <= (
                1.5*(coarse['omitted_mass_union_bound']+fine['omitted_mass_union_bound'])+
                coarse['numerical_budget']+fine['numerical_budget'])
            assert fine['gene']['TV_lower'] == 0
            if label == 'heterogeneous':
                assert fine['bin']['TV_lower'] > 0
                assert np.max(np.abs(np.array(fine['joint_low'])-np.array([[float(v) for v in row] for row in joint]))) <= fine['numerical_budget']
            else:
                assert fine['bin']['TV_lower'] == 0
            examples.append(dict(label=label, means=means.tolist(), sizes=sizes.tolist(), split_size=split_size,
                coarse=coarse, refined=fine))
        contexts = []
        for seed in (20261108, 20261109, 20261110):
            values = np.load(inputs[f'{seed}:parameters'])
            mapping = pd.read_csv(inputs[f'{seed}:mapping'], sep='\t')
            selected_bins = sorted(set(mapping.loc[mapping.chromosome.eq('chr12'), 'bin']))[:2]
            assert len(selected_bins) == 2
            exposure = float(np.median(values['L']))
            phi = float(values['phi'])
            assert phi == .1
            for k in range(3):
                rate = values['class_rates'][k]*(exposure/float(values['depth']))**values['gamma']
                means = exposure*rate/rate.sum()
                for j in selected_bins:
                    selected = values['source_columns'][values['gene_bins'] == j]
                    mu, size = means[selected], np.full(len(selected), 1/phi)
                    moment_size = mu.sum()**2/np.sum(mu**2/size)
                    q = mu/(size+mu)
                    s1, s2 = float(np.sum(size*q)), float(np.sum(size*q*q))
                    f0 = float(np.exp(-np.sum(size*np.log1p(mu/size))))
                    f1, f2 = f0*s1, f0*(s1*s1+s2)/2
                    low = np.array([[f0, f1/2], [f1/2, f2*moment_size/(2*(moment_size+1))]])
                    independent = distribution(mu, size, 2)
                    assert np.max(np.abs(independent-[f0, f1, f2])) <= 128*len(mu)*np.finfo(float).eps
                    odds = float(low[0, 0]*low[1, 1]/(low[0, 1]*low[1, 0]))
                    direct = moment_size/(moment_size+1)*(1+s2/(s1*s1))
                    assert abs(odds-direct) <= 128*len(mu)*np.finfo(float).eps
                    contexts.append(dict(seed=seed, class_index=k, bin=int(j), exposure=exposure,
                        genes=len(mu), conditional_mean=float(mu.sum()), moment_matched_size=float(moment_size),
                        low_joint=low.tolist(), low_joint_odds_ratio=odds,
                        scope='Low-count conditional diagnostic, not total variation or specimen error'))
        (a.out/'report.json').write_text(json.dumps(dict(all_checks_passed=True,
            exact_counterexample=dict(mean=str(M), variance=str(V), split_size=str(R), covariance=str(covariance),
                low_joint=[[str(v) for v in row] for row in joint], determinant=str(determinant)),
            examples=examples, original_contexts=contexts, biological_specimens=0,
            production_changed=False, scientific_criteria_changed=False,
            scope='Fixed known noise; zero covariance does not establish independence; gene-native independence is conditional'),
            indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU exact Fraction, finite-grid tail/refinement and independent quadrature')


if __name__ == '__main__':
    main()
