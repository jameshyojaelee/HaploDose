"""Independent count/variance/likelihood references for gene event proposals."""
import argparse
from fractions import Fraction
import importlib.util
from pathlib import Path

import numpy as np
from scipy.special import expit, logsumexp
from scipy.stats import nbinom

from haplodose.clones import _event_mixture_objective
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    path = Path(__file__).with_name('gene_event_proposal.py')
    spec = importlib.util.spec_from_file_location('checked_gene_events',path)
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    run = start_run(args.out,dict(script=Path(__file__),model=path,design=args.task/'DESIGN.md'),{},
        dict(scope='Numerical invariants, not independent biological validation'),mode='gene_event_numerics')
    try:
        counts = np.array([[1.,4.,2.,0.],[3.,8.,4.,0.]])
        gbin = np.array([0,0,1,1])
        sizes = np.full(4,2.)
        fractions, bin_sizes, supported = model.composition(counts,gbin,sizes,2)
        np.testing.assert_array_equal(fractions,[float(Fraction(1,4)),float(Fraction(3,4)),1.,0.])
        np.testing.assert_array_equal(bin_sizes,[float(Fraction(16,5)),2.])
        assert supported.all()
        expected_bins = np.array([[6.,8.],[3.,4.]])  # Order: bin1, then bin0.
        O,E,r,genes,local = model.region_arrays(counts,gbin,fractions,sizes,[1,0],expected_bins)
        np.testing.assert_array_equal(genes,[0,1,2])
        np.testing.assert_array_equal(local,[1,1,0])
        for j,b in enumerate([1,0]):
            np.testing.assert_array_equal(O[:,local==j].sum(axis=1),counts[:,gbin==b].sum(axis=1))
            np.testing.assert_array_equal(E[:,local==j].sum(axis=1),expected_bins[:,j])
            variance = (E[:,local==j]+E[:,local==j]**2/r[local==j]).sum(axis=1)
            np.testing.assert_array_equal(variance,expected_bins[:,j]+expected_bins[:,j]**2/bin_sizes[b])
        x = np.array([-.2,.3,-.4])
        null_eta = .1
        value,gradient,q = _event_mixture_objective(x,O,E,r,null_eta)
        logweights = np.log([1-expit(x[2]),expit(x[2])])
        scores = np.array([nbinom.logpmf(O,r,r/(r+E*np.exp(eta))).sum(axis=1) for eta in x[:2]]).T+logweights
        native_null = nbinom.logpmf(O,r,r/(r+E*np.exp(null_eta))).sum(axis=1)
        direct = -float(np.sum(logsumexp(scores,axis=1)-native_null))
        direct_q = np.exp(scores[:,1]-logsumexp(scores,axis=1))
        budget = (128*O.size+32)*np.finfo(float).eps*(1+abs(direct)+np.abs(scores).sum())
        assert abs(value-direct) <= budget
        assert np.max(np.abs(q-direct_q)) <= budget
        permutation = np.array([2,0,1])
        other = _event_mixture_objective(x,O[:,permutation],E[:,permutation],r[permutation],null_eta)
        assert abs(other[0]-value) <= budget and np.max(np.abs(other[2]-q)) <= budget
        # Two different gene allocations have the same bin sum. Exact rational
        # geometric PMFs independently show their class likelihood ratios differ.
        allocation = np.array([[0.,1.],[1.,0.]])
        mu = [Fraction(1),Fraction(4)]
        ratios=[]
        for row in allocation.astype(int):
            probabilities=[]
            for theta in (Fraction(1),Fraction(3,2)):
                p=Fraction(1)
                for y,m in zip(row,mu):
                    scaled=theta*m
                    p *= Fraction(1,1)/(1+scaled)*(scaled/(1+scaled))**int(y)
                probabilities.append(p)
            ratios.append(probabilities[1]/probabilities[0])
        _,_,gene_q=_event_mixture_objective(np.array([0.,np.log(1.5),0.]),allocation,
                                          np.array([[1.,4.],[1.,4.]]),np.ones(2),0.)
        for i,ratio in enumerate(ratios):
            assert abs(gene_q[i]-float(ratio/(1+ratio))) <= budget
        assert ratios[0]!=ratios[1] and gene_q[0]!=gene_q[1]
        # An unobserved bin remains unresolved instead of receiving a floor.
        _,empty_sizes,empty_support=model.composition([[1.,0.]],np.array([0,1]),[2.,2.],2)
        assert not empty_support[1] and np.isnan(empty_sizes[1])
        write_json(dict(all_checks_passed=True,count_expectation_variance_conserved=True,
            independent_scipy_likelihood_discrepancy=abs(value-direct),operation_rounding_budget=float(budget),
            gene_permutation_invariant=True,same_sum_distinct_rational_likelihood_ratios=[str(r) for r in ratios],
            missing_bin_abstention=True),args.out/'report.json')
        finish_run(args.out,run,backend='CPU Fraction moments and independent SciPy NB likelihood')
    except BaseException as error:
        finish_run(args.out,run,error=error)
        raise


if __name__=='__main__':
    main()
