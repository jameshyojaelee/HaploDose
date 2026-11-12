"""Correlation of exact arm summaries of finite stored RNA/DNA values.

The observation is the supplied binary floating-point value. Mean/median
aggregation and average-tie ranks use rational arithmetic so equal summaries
remain tied regardless of arm length, input dtype or summation order. Display
correlations use 80 decimal digits; no numerical tolerance defines a win.
"""
from decimal import Decimal, localcontext
from fractions import Fraction
from math import isfinite
from statistics import median


def _ranks(values):
    ordered = sorted(set(values))
    counts = {value: values.count(value) for value in ordered}
    mapping, start = {}, 1
    for value in ordered:
        count = counts[value]
        mapping[value] = Fraction(2 * start + count - 1, 2)
        start += count
    return [mapping[value] for value in values]


def _correlation(x, y):
    n = len(x)
    covariance = n * sum(a * b for a, b in zip(x, y)) - sum(x) * sum(y)
    vx = n * sum(a * a for a in x) - sum(x) ** 2
    vy = n * sum(b * b for b in y) - sum(y) ** 2
    if vy == 0:
        return float('nan')
    if vx == 0:
        return 0.
    with localcontext() as context:
        context.prec = 80
        def decimal(value):
            return Decimal(value.numerator) / Decimal(value.denominator)
        return float(decimal(covariance) / (decimal(vx) * decimal(vy)).sqrt())


def aggregate_metrics(values, truth, groups, aggregation):
    """Score exact unweighted arm means/medians on already fixed support."""
    if aggregation not in ('mean', 'median'):
        raise ValueError('Expected mean or median aggregation')
    if not (len(values) == len(truth) == len(groups)):
        raise ValueError('Aligned values, truth and arm labels required')
    grouped = {}
    for value, target, group in zip(values, truth, groups):
        if not isfinite(float(value)) or not isfinite(float(target)):
            raise ValueError('Finite values on fixed callable support required')
        if group is None or not str(group).strip():
            raise ValueError('Nonempty arm labels required')
        x, y = grouped.setdefault(str(group), ([], []))
        x.append(Fraction(float(value)))
        y.append(Fraction(float(target)))
    x, y = [], []
    for group in sorted(grouped):
        first, second = grouped[group]
        reduce = (lambda v: sum(v) / len(v)) if aggregation == 'mean' else median
        x.append(reduce(first))
        y.append(reduce(second))
    if len(x) < 4:
        return dict(n=len(x), pearson=float('nan'), spearman=float('nan'))
    return dict(n=len(x), pearson=_correlation(x, y), spearman=_correlation(_ranks(x), _ranks(y)))
