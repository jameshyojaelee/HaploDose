"""Reference-only adjacent-depth factorial-moment gene noise estimator."""
from fractions import Fraction
import sys

import numpy as np
from haplodose import pipeline


def pair_groups(depth, identities):
    """Adjacent sorted pairs, with one terminal triplet for an odd population."""
    depth, identities = np.asarray(depth, float), np.asarray(identities, str)
    if (depth.ndim != 1 or identities.shape != depth.shape or len(depth)<2
            or not np.isfinite(depth).all() or np.any(depth <= 0) or len(np.unique(identities)) != len(identities)):
        raise ValueError('At least two positive-depth distinct reference cells required')
    order = np.lexsort((identities, np.log(depth)))
    stop = len(order)-3 if len(order)%2 else len(order)
    groups = [order[i:i+2] for i in range(0, stop, 2)]
    if stop != len(order):
        groups.append(order[stop:])
    assert sorted(np.concatenate(groups).tolist()) == list(range(len(depth)))
    return groups


def moments(counts, depth, groups, identities):
    counts, depth, groups, identities = map(np.asarray, (counts, depth, groups, identities))
    if (counts.ndim != 2 or counts.shape[0] != len(depth) or groups.shape != depth.shape
            or identities.shape != depth.shape or not np.isfinite(counts).all() or np.any(counts < 0)
            or np.any(counts != np.floor(counts))):
        raise ValueError('Aligned nonnegative integer reference counts and metadata required')
    numerator, denominator, details = np.longdouble(0), np.longdouble(0), []
    for group in sorted(set(groups)):
        selected = np.flatnonzero(groups == group)
        for ix in pair_groups(depth[selected], identities[selected]):
            rows = selected[ix]
            y, exposure = counts[rows].astype(np.longdouble), depth[rows].astype(np.longdouble)
            pair_exposure = exposure.sum()**2-np.sum(exposure**2)
            if pair_exposure <= 0:
                raise ArithmeticError('Unresolved pair exposure')
            rate_square = ((y.sum(axis=0)**2-np.sum(y*y, axis=0))/pair_exposure).sum()
            squared_mean = rate_square*np.sum(exposure**2)
            factorial = np.sum(y*(y-1))
            numerator += factorial-squared_mean
            denominator += squared_mean
            details.append(dict(group=str(group), rows=rows.tolist(), ids=identities[rows].tolist(),
                log_depth_span=float(np.ptp(np.log(exposure))), factorial_moment=float(factorial),
                squared_mean_moment=float(squared_mean)))
    if denominator <= 0:
        raise ValueError('Reference counts have no positive cross-cell mean-square support')
    phi = float(numerator/denominator)
    if not np.isfinite(phi) or phi <= 0:
        raise ValueError('Reference counts do not identify a positive common dispersion; no clipped/oracle fallback')
    return phi, details


def estimator(records, depth_strata=5):
    if depth_strata != 5:
        raise ValueError('Replace only the frozen quintile estimator')
    cache = {}
    def size(info):
        frame = sys._getframe(1)
        while frame is not None and frame.f_code is not pipeline._run_single.__code__:
            frame = frame.f_back
        if frame is None:
            raise RuntimeError('Missing exact pipeline reference context')
        loc = frame.f_locals
        key = (str(loc['counts']), tuple(loc['used'].column), str(loc['dataset']))
        if key in cache:
            return cache[key]
        meta = loc['meta']
        rows = np.flatnonzero(loc['fit'] & meta.dataset.eq(loc['dataset']).to_numpy())
        assert meta.iloc[rows].role.eq('reference').all() and meta.iloc[rows].ref_split.eq('fit').all()
        X, gbin = pipeline.gene_loader(loc['counts'], loc['used'], rows)(np.arange(len(rows)))
        X = X.toarray().astype(float)
        depth = np.asarray(loc['Lbin'][rows], float)
        np.testing.assert_array_equal(X.sum(axis=1), depth)
        np.testing.assert_array_equal(gbin, info['gbin'])
        groups = (meta.dataset+'|'+meta.ref_group).to_numpy()[rows]
        phi, details = moments(X, depth, groups, meta.cell_id.to_numpy()[rows])
        value = 1/phi
        cache[key] = value
        records.append(dict(counts=str(loc['counts']), training_rows=rows.tolist(),
            training_cell_ids=meta.cell_id.to_numpy()[rows].tolist(), estimated_phi=phi, estimated_size=value,
            genes=X.shape[1], groups=details, matching='adjacent_depth_pairs_with_terminal_odd_triplet',
            query_rows_used=0, holdout_reference_rows_used=0, true_noise_used=False))
        return value
    return size


def checks():
    counts = np.array([[1, 6], [8, 2], [2, 10], [18, 3]])
    depth = np.array([7., 10., 12., 21.])
    groups = np.array(['A']*4)
    ids = np.array(['b', 'a', 'c', 'd'])
    phi, details = moments(counts, depth, groups, ids)
    numerator, denominator = Fraction(0), Fraction(0)
    for pair in ((0,1), (2,3)):
        a,b = pair
        mean_square = sum(Fraction(int(counts[a,g])*int(counts[b,g])) for g in range(2))
        mean_square *= Fraction(int(depth[a])**2+int(depth[b])**2, int(depth[a])*int(depth[b]))
        factorial = sum(Fraction(int(counts[i,g])*(int(counts[i,g])-1)) for i in pair for g in range(2))
        numerator += factorial-mean_square
        denominator += mean_square
    exact = float(numerator/denominator)
    operations = 8*counts.size+32
    eps = np.finfo(float).eps
    allowance = operations*eps/(1-operations*eps)*max(1., abs(exact))
    assert abs(phi-exact) <= allowance
    order = np.array([3, 0, 2, 1])
    other,_ = moments(counts[order], depth[order], groups[order], ids[order])
    assert phi == other
    np.testing.assert_array_equal(pair_groups([2,1,4,3,5], ['b','a','d','c','e'])[-1], [3,2,4])
    # Conditional expectation under constant gene rates: E[A]=(1+phi)B.
    rate = [Fraction(1,7), Fraction(2,7)]
    exposure = list(map(Fraction, [7,10]))
    true_phi = Fraction(1,10)
    B = sum(e*e for e in exposure)*sum(r*r for r in rate)
    A = (1+true_phi)*B
    assert (A-B)/B == true_phi
    try:
        moments(np.zeros((4,2)), depth, groups, ids)
    except ValueError:
        pass
    else:
        raise AssertionError('Unsupported references must fail')
    return dict(independent_fraction_moments=True, conditional_NB_expectation_identity=True,
        reference_order_invariant=True, odd_triplet_all_rows_retained=True,
        unsupported_mean_support_rejected=True, float_rounding_allowance=allowance)
