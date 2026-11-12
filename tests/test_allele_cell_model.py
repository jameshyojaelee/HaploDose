"""Independent allele observation-model contrasts through selection paths."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import sparse, stats
from scipy.special import expit, gammaln

from haplodose.alleles import (AlleleParams, SpecimenAlleles, balance_loglik_ratio,
                              balanced_cells, cell_loglik)
from haplodose.clones import CloneParams, SpecimenModel


def test_balance_contrast_matches_independent_beta_binomial():
    a = sparse.csr_matrix([[19.,800.],[1.,1400.],[0.,0.]])
    d = sparse.csr_matrix([[20.,2000.],[20.,2000.],[0.,0.]])
    b = np.array([.77,.77])
    kappa = 2.
    expected = (stats.betabinom.logpmf(a.toarray(),d.toarray(),kappa*b,kappa*(1-b))
                -stats.betabinom.logpmf(a.toarray(),d.toarray(),kappa*.5,kappa*.5)).sum(axis=1)
    got = balance_loglik_ratio(a,d,b,kappa)
    # PMF subtraction cancels count-dependent log-gamma terms. Scale the
    # arithmetic allowance to those terms, rather than using a universal epsilon.
    allowance = 64*np.finfo(float).eps*np.sum(1+abs(gammaln(d.toarray()+1)))
    np.testing.assert_allclose(got,expected,atol=allowance,rtol=0)
    assert got.dtype == np.float64 and got[2] == 0.
    flipped = balance_loglik_ratio(d-a,d,1-b,kappa)
    np.testing.assert_allclose(flipped,got,atol=allowance,rtol=0)
    # More reads at the same fraction do not remove finite-concentration ASE.
    deep = balance_loglik_ratio(10*a,10*d,b,kappa)
    assert abs(deep[0]) < 1.


def test_legacy_contrast_and_membership_precision_are_unchanged():
    a = sparse.csr_matrix([[19.,800.],[1.,1400.]])
    d = sparse.csr_matrix([[20.,2000.],[20.,2000.]])
    b = np.array([.77,.77])
    expected = np.asarray(a@np.log(b/.5)+(d-a)@np.log((1-b)/.5)).ravel()
    np.testing.assert_array_equal(balance_loglik_ratio(a,d,b),expected)
    assert cell_loglik(a,d,b[None],kappa=2.).dtype == np.float32


@pytest.mark.parametrize('kappa',[0.,-1.,np.inf,np.nan])
def test_invalid_concentration_cannot_enter_selection(kappa):
    with pytest.raises(ValueError,match='finite and positive'):
        balance_loglik_ratio(sparse.csr_matrix([[1.]]),sparse.csr_matrix([[2.]]),np.array([.77]),kappa)


@pytest.mark.parametrize('kappa',[None,2.])
def test_allele_proposals_use_selected_cell_model(monkeypatch,kappa):
    a = np.array([[19,0,800,0],[1,0,1400,0],[10,0,1000,0]],float)
    d = np.array([[20,0,2000,0]]*3,float)
    snps = pd.DataFrame(dict(chromosome='chr1',pos=[1,2,3,4],cM=[0.,.01,.02,.03]))
    alleles = SpecimenAlleles(a,d,snps,50.,AlleleParams(cell_kappa=kappa))
    post = np.zeros((4,9));post[:,5] = 1.  # Oriented b=.8, then existing .9 shrink -> .77.
    monkeypatch.setattr(alleles,'runs',lambda q,n:[dict(resolved=True,first_snp=0,last_snp=3,
                                                       chromosome='chr1',start=0,end=4)])
    model = SimpleNamespace(alleles=alleles,allele_post=post[None],Dk=d.sum(axis=0)[None],
                            params=CloneParams(min_cells=1))
    proposals = SpecimenModel.allele_splits(model,0,np.arange(3))
    if kappa is None:
        contrast = (stats.binom.logpmf(a[:,::2],d[:,::2],.77)
                    -stats.binom.logpmf(a[:,::2],d[:,::2],.5)).sum(axis=1)
    else:
        contrast = (stats.betabinom.logpmf(a[:,::2],d[:,::2],2*.77,2*.23)
                    -stats.betabinom.logpmf(a[:,::2],d[:,::2],1.,1.)).sum(axis=1)
    assert len(proposals) == 1
    np.testing.assert_array_equal(proposals[0][0],contrast>0)
    np.testing.assert_array_equal(proposals[0][0],[0,1,0] if kappa is None else [1,0,0])


@pytest.mark.parametrize('kappa',[None,2.])
def test_balanced_cell_selection_respects_overdispersion(monkeypatch,kappa):
    snps = pd.DataFrame(dict(chromosome=['chr1'],pos=[1],cM=[0.]))
    alleles = SpecimenAlleles([[500.]],[[1000.]],snps,50.,AlleleParams(cell_kappa=kappa))
    post = np.zeros((1,1,9));post[:,:,5] = 1.
    monkeypatch.setattr(alleles,'fit',lambda a,d:(post,np.array([[.8]])))
    monkeypatch.setattr(alleles,'runs',lambda q,n:[dict(resolved=True,theta=.3,first_snp=0,last_snp=0)])
    if kappa is None:
        llr = stats.binom.logpmf(500,1000,.77)-stats.binom.logpmf(500,1000,.5)
    else:
        llr = stats.betabinom.logpmf(500,1000,2*.77,2*.23)-stats.betabinom.logpmf(500,1000,1.,1.)
    expected = [0] if expit(-llr)>=.99 else []
    assert balanced_cells(alleles)[0].tolist() == expected
    assert expected == ([0] if kappa is None else [])
