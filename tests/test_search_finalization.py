"""A last-cycle split must face the same merge rule as earlier splits."""
from types import SimpleNamespace
import numpy as np
import pytest
from haplodose.clones import _finalize_search


def state_machine(target):
    """Hand-specified unsupported states: remove one group until target."""
    calls=[]
    def em(doses,pi):
        n=len(doses);calls.append(('refit',n))
        # Unique per-state arrays expose accidental return of stale pools or R.
        return doses,np.full(n,1/n),np.full((5,n),1/n),np.full((5,n),n),np.full((n,2),n),np.full(n,n),np.full((n,2),-n)
    def merge(R,doses,pi):
        n=len(doses);calls.append(('check',n))
        assert R.shape==(5,n)
        if n>target:return doses[:-1],np.full(n-1,1/(n-1)),True
        return doses,pi,False
    return SimpleNamespace(em=em,merge_unsupported=merge,log=[]),calls


@pytest.mark.parametrize('initial,expected',[
    (2,[('refit',2),('check',2)]),
    (3,[('refit',3),('check',3),('refit',2),('check',2)]),
    (4,[('refit',4),('check',4),('refit',3),('check',3),('refit',2),('check',2)]),
])
def test_unsupported_terminal_groups_are_refitted_and_rechecked(initial,expected):
    model,calls=state_machine(2)
    doses,pi,R,LL,Y,L,XL=_finalize_search(model,[np.full(2,k) for k in range(initial)],np.full(initial,1/initial))
    assert calls==expected
    assert len(doses)==2
    np.testing.assert_array_equal(pi,[.5,.5])
    np.testing.assert_array_equal(R,np.full((5,2),.5))
    np.testing.assert_array_equal(LL,np.full((5,2),2))
    np.testing.assert_array_equal(Y,np.full((2,2),2))
    np.testing.assert_array_equal(L,[2,2])
    np.testing.assert_array_equal(XL,np.full((2,2),-2))
    assert model.log[-1]['terminal_merge_rounds']==initial-2
    assert model.log[-1]['reason']=='cycle_budget_exhausted'


def test_nondecreasing_merge_cannot_create_an_unbounded_finalization_loop():
    model,_=state_machine(2)
    model.merge_unsupported=lambda R,doses,pi:(doses,pi,True)
    with pytest.raises(RuntimeError,match='reduce'):
        _finalize_search(model,[np.zeros(2),np.ones(2)],np.array([.5,.5]))
