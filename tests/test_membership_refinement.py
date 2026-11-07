"""Stable hard assignments are not a convergence certificate for soft refinement."""
from types import SimpleNamespace
import numpy as np
from scipy.special import softmax
from haplodose.clones import SpecimenModel


def test_stable_map_labels_do_not_stop_changing_profiles():
    # Isolate the refinement controller with an independently solvable contraction:
    # d[t+1]=(d[t]+1)/2, hence d[t]=1-2**(-t). Every MAP label remains tumour,
    # although its probability and the profile move substantially after round 2.
    m = SpecimenModel.__new__(SpecimenModel)
    m.params = SimpleNamespace(outer_rounds=8, normal_prior=.05, min_cells=0)
    m.J = 1
    m.alleles = None
    m.R_prev = None
    m.fixed_gamma = True
    m.log = []
    m.lam = 0.
    seen = []

    def loglik(doses):
        m.current = float(doses[1][0])
        return np.tile([0., 1 + m.current], (10, 1))

    def profile(*args, **kwargs):
        value = (m.current + 1) / 2
        seen.append(value)
        return dict(dose=np.array([[0.], [value]]), E=np.ones((2,1)),
                    post=np.ones((2,1,1)), path=np.zeros((2,1)), tau=np.ones(2))

    m.loglik = loglik
    m.posterior = lambda LL, pi: softmax(LL + np.log(pi)[None], axis=1)
    m.pooled = lambda R: (R.sum(axis=0)[:,None], R.sum(axis=0), np.zeros((2,2)))
    m.shared_profiles = profile
    m.has_events = lambda *args: True
    m.allele_events = lambda k: False
    m._keep_alleles = lambda keep: None
    m._update_offset = lambda *args: None
    m._update_phi = lambda *args: None
    m.update_scores = lambda *args: None
    doses, _, R, *_ = m.em([np.zeros(1),np.zeros(1)],np.array([.5,.5]), stop_on_map=False)
    assert doses[1][0] == 1 - 2**(-m.params.outer_rounds)
    assert len(seen) == m.params.outer_rounds
    assert np.all(R.argmax(axis=1) == 1)
    assert m.log[-1]['convergence_status'] == 'not_certified_budgeted'
    # Preserve the historical default as an explicitly labelled heuristic.
    # Its second iterate is not the eight-round result or the fixed point.
    doses, *_ = m.em([np.zeros(1),np.zeros(1)],np.array([.5,.5]))
    assert doses[1][0] == .75
    assert m.log[-1]['convergence_status'] == 'not_certified_map_stability'
    assert m.log[-1]['refinement_rounds'] == 2


def test_merge_reports_assigned_pool_uncertainty_without_equivalence_claim():
    m=SpecimenModel.__new__(SpecimenModel)
    m.params=SimpleNamespace(min_cells=30)
    m.J=1;m.thin=False;m.Ak=None;m.allele_post=None;m.graph=None;m.log=[]
    R=np.tile([.1,.46,.44],(60,1));R[30:,1:]=R[30:,1:][:,::-1]
    m.matched_differences=lambda labels,pairs:[dict(segments=[],distance=0.,testable=True)]
    m.supported=lambda *args,**kwargs:[]
    m.pooled=lambda R:(np.ones((2,1)),np.ones(2),np.zeros((2,2)))
    m.profiles=lambda *args:dict(dose=np.zeros((2,1)))
    _,_,merged=m.merge_unsupported(R,[np.zeros(1),np.ones(1),np.full(1,2.)],np.array([.1,.45,.45]))
    assert merged
    event=m.log[0]
    q=.46/.9
    expected=-q*np.log(q)-(1-q)*np.log1p(-q)
    assert abs(event['mean_pair_entropy_nats']-expected)<=16*np.finfo(float).eps
    assert event['criterion_scope']=='assigned_cell_pools'
    assert event['interpretation']=='no_resolved_difference_not_genomic_equivalence'
