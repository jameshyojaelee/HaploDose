"""Research-only profile-baseline conditioning, restored outside each profile call."""
from contextlib import contextmanager

import numpy as np
from scipy.special import softmax
from haplodose.clones import SpecimenModel


@contextmanager
def conditioned_profiles(trace, contexts):
    original_profile, original_expected = SpecimenModel.profiles, SpecimenModel.expected
    active, identities = {}, {}

    def identity(model):
        if id(model) not in identities:
            identities[id(model)] = (len(identities), model)
        return identities[id(model)][0]

    def expected(model, Lk, XLk):
        result = original_expected(model, Lk, XLk)
        if id(model) in active:
            base, offset, gamma = active[id(model)]
            key = (identity(model), bool(np.any(offset != offset[0])))
            if key not in contexts:
                contexts[key] = dict(base=base.copy(), offset=offset.copy(), gamma=model.gamma.copy(),
                    original_gamma=gamma.copy(),
                    Lk=np.atleast_1d(Lk).copy(), XLk=np.atleast_2d(XLk).copy(),
                    profile_base=model.base.copy(), expected=result.copy())
        return result

    def profiles(model, *args, **kwargs):
        if id(model) in active:
            raise RuntimeError('Nested profile conditioning is unsupported')
        base, offset, gamma = model.base, model.offset.copy(), model.gamma
        before = base.tobytes()
        before_gamma = gamma.tobytes()
        active[id(model)] = (base, offset, gamma)
        if np.any(offset != offset[0]):
            model.base = softmax(np.log(base)+offset)
            model.gamma = gamma-float(model.base@gamma)
        try:
            return original_profile(model, *args, **kwargs)
        finally:
            model.base = base
            model.gamma = gamma
            del active[id(model)]
            assert model.base.tobytes() == before
            assert model.gamma.tobytes() == before_gamma
            assert np.array_equal(model.offset, offset)
            trace.append(dict(model=identity(model), bins=len(base),
                nonconstant_offset=bool(np.any(offset != offset[0])),
                original_base_restored=True, original_gamma_restored=True, offset_unchanged=True))

    SpecimenModel.profiles, SpecimenModel.expected = profiles, expected
    try:
        yield
    finally:
        SpecimenModel.profiles, SpecimenModel.expected = original_profile, original_expected
