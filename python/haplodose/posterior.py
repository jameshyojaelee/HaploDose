"""Descriptive RNA profiles averaged over working regional boundaries.

The Gaussian weights and state grid define a smoother, not calibrated DNA
states or interval estimates. Existing clone/bulk inference defaults do not
call this module. Prepared-data callers preserve gene/reference correspondence.
"""
from dataclasses import asdict
import math

import numpy as np
import pandas as pd
from scipy.special import logsumexp, ndtri

from . import bulk, profile


def _transition(prior, change):
    prior=np.asarray(prior,float)
    if np.any(prior<=0) or not np.isclose(prior.sum(),1.) or not 0<change<1:
        raise ValueError('Invalid state prior/transition probability')
    prior=prior/prior.sum()
    T=change*prior[None,:]/(1-prior[:,None])
    np.fill_diagonal(T,1-change)
    return T


def _chain(emission, logT, initial, terminal):
    """Double-precision forward/backward and MAP path for one interval."""
    e=np.asarray(emission,float);n,k=e.shape
    forward=np.empty_like(e);forward[0]=initial+e[0]
    best=forward[0].copy();back=np.zeros((n,k),int)
    for j in range(1,n):
        forward[j]=e[j]+logsumexp(forward[j-1][:,None]+logT,axis=0)
        scores=best[:,None]+logT
        back[j]=scores.argmax(axis=0);best=e[j]+scores.max(axis=0)
    evidence=float(logsumexp(forward[-1]+terminal))
    backward=np.empty_like(e);backward[-1]=terminal
    for j in range(n-2,-1,-1):
        backward[j]=logsumexp(logT+e[j+1][None,:]+backward[j+1][None,:],axis=1)
    posterior=np.exp(forward+backward-evidence)
    path=np.empty(n,int);path[-1]=np.argmax(best+terminal)
    for j in range(n-2,-1,-1):path[j]=back[j+1,path[j+1]]
    return posterior,evidence,path


def _weighted_median(values, weights):
    """Exact dyadic weight comparisons, including half-mass ties."""
    values, weights = np.asarray(values, float), np.asarray(weights, float)
    if (values.ndim != 1 or weights.shape != values.shape or not len(values)
            or not np.isfinite(values).all() or not np.isfinite(weights).all() or np.any(weights <= 0)):
        raise ValueError('Finite effects and strictly positive aligned weights required')
    order = np.argsort(values, kind='stable')
    values, weights = values[order], weights[order]
    dyadic = [float(w).as_integer_ratio() for w in weights]
    exponent = max(d.bit_length()-1 for _, d in dyadic)
    integers = [n << (exponent-(d.bit_length()-1)) for n, d in dyadic]
    total, cumulative = sum(integers), 0
    for j, weight in enumerate(integers):
        cumulative += weight
        if 2*cumulative >= total:
            return float((values[j]+values[j+1])/2) if 2*cumulative == total else float(values[j])
    raise ArithmeticError('Weighted median mass was not conserved')


def _gene_variance(spread, expected, weighting):
    """Positive descriptive weighting scale, without an inferential SE claim."""
    spread, expected = np.asarray(spread, float), np.asarray(expected, float)
    if (spread.shape != expected.shape or not np.isfinite(spread).all() or np.any(spread < 0)
            or not np.isfinite(expected).all() or np.any(expected < 0)
            or weighting not in ('reference', 'count')):
        raise ValueError('Aligned finite reference spread and supported expected counts required')
    return (spread if weighting == 'reference' else 0)+1/(expected+.5)


def _normal_centered_effect(values, normal_effects, gene_bins):
    """Relative RNA against typical normals, with the existing bin-median centre."""
    values, normal_effects, gene_bins = map(np.asarray, (values, normal_effects, gene_bins))
    if (values.ndim != 1 or normal_effects.ndim != 2 or normal_effects.shape[1] != len(values)
            or not len(normal_effects) or gene_bins.shape != values.shape or gene_bins.dtype.kind not in 'iu'
            or not np.array_equal(np.unique(gene_bins), np.arange(int(gene_bins.max())+1))
            or not np.isfinite(values).all() or not np.isfinite(normal_effects).all()):
        raise ValueError('Finite aligned query/normal gene effects and complete integer bin IDs required')
    location = np.median(normal_effects, axis=0)
    effect = values-location
    centre = float(np.median([np.median(effect[gene_bins == j]) for j in np.unique(gene_bins)]))
    return effect-centre, location, centre


def _bin_observations(effects,expected,variance,gene_bins):
    y,e,v,bins=map(np.asarray,(effects,expected,variance,gene_bins))
    if (y.ndim!=1 or not len(y) or any(x.shape!=y.shape for x in (e,v,bins))
            or bins.dtype.kind not in 'iu' or not np.array_equal(np.unique(bins),np.arange(int(bins.max())+1))
            or not all(np.isfinite(x).all() for x in (y,e,v)) or np.any(e<0) or np.any(v<=0)):
        raise ValueError('Finite aligned effects, expectations, positive scales and complete bin IDs required')
    median=_weighted_median
    values=[];sampling=[]
    for j in range(int(bins.max())+1):
        keep=bins==j;a=1/np.sqrt(v[keep]);density=a*np.sqrt(e[keep]+.5)
        values.append(median(y[keep],a))
        sampling.append((np.pi/2)*math.fsum(map(float,a*a))/math.fsum(map(float,density))**2)
    return np.asarray(values),np.asarray(sampling)


def _train(normal_effects,normal_expected,gene_bins):
    effects,expected=np.asarray(normal_effects,float),np.asarray(normal_expected,float)
    if effects.ndim!=2 or expected.shape!=effects.shape or len(effects)<3:
        raise ValueError('At least three aligned normal profiles required')
    rows=[]
    for i in range(len(effects)):
        other=np.delete(effects,i,axis=0)
        adjusted,centre,_=_normal_centered_effect(effects[i],other,gene_bins)
        spread=(np.median(abs(other-centre),axis=0)/ndtri(.75))**2
        variance=_gene_variance(spread,expected[i],'reference')
        rows.append(_bin_observations(adjusted,expected[i],variance,gene_bins)[0])
    rows=np.asarray(rows);location=np.median(rows,axis=0)
    variance=(np.median(abs(rows-location),axis=0)/ndtri(.75))**2
    return dict(normal_bins=rows,reference_variance=variance,
        status='Conditional working spread; no calibrated sampling or DNA-state probabilities')


def _chain_blocks(emission,reset,prior,change):
    emission,reset,prior=np.asarray(emission,float),np.asarray(reset,bool),np.asarray(prior,float)
    if emission.ndim!=2 or reset.shape!=(len(emission),) or not len(reset) or not reset[0]:
        raise ValueError('Aligned emissions and genomic resets required')
    prior=prior/prior.sum()
    logT=np.log(_transition(prior,change));initial=np.log(prior)
    starts=np.flatnonzero(reset);stops=np.r_[starts[1:],len(reset)]
    post=np.empty_like(emission);path=np.empty(len(reset),int);scores=[]
    for a,b in zip(starts,stops):
        post[a:b],value,path[a:b]=_chain(emission[a:b],logT,initial,np.zeros(len(prior)))
        scores.append(value)
    return post,math.fsum(scores),path


def _fit(effects,expected,gene_variance,gene_bins,reference_variance,reset,params,*,centre_mode='zero'):
    if centre_mode not in ('zero','bin-median'):raise ValueError('Unsupported boundary intercept')
    values,sampling=_bin_observations(effects,expected,gene_variance,gene_bins)
    variance=np.asarray(reference_variance,float)+sampling;reset=np.asarray(reset,bool)
    if (variance.shape!=values.shape or reset.shape!=values.shape or not reset[0]
            or not np.isfinite(variance).all() or np.any(variance<=0)):
        raise ValueError('Finite positive bin scales and complete genomic reset correspondence required')
    inside=~reset[1:];difference=np.diff(values)[inside]
    if len(difference)<10:tau=params.tau_floor
    else:
        mad=np.median(abs(difference-np.median(difference)))/ndtri(.75)
        tau=math.sqrt(max(params.tau_floor**2,(mad**2-float(np.median((variance[:-1]+variance[1:])[inside])))/2))
    variance=variance+tau**2;fits=[]
    centre=0. if centre_mode=='zero' else float(np.median(values))
    for beta in bulk.BETA_GRID:
        if not params.beta_bounds[0]<=beta<=params.beta_bounds[1]:continue
        means=beta*np.asarray(profile.LOG_DOSE,float)
        emission=-.5*(np.log(2*np.pi*variance)[:,None]+(values[:,None]-centre-means)**2/variance[:,None])
        post,score,path=_chain_blocks(emission,reset,profile.STATE_PRIOR,params.p_change)
        fits.append(dict(beta=beta,score=score,path=path,posterior=post))
    if not fits:raise ValueError('No existing bulk beta grid value satisfies the declared bounds')
    best=max(range(len(fits)),key=lambda j:fits[j]['score']);selected=fits[best]
    operations=128*len(values)*len(profile.CN)+256
    allowance=operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)*max(1.,max(abs(x['score']) for x in fits))
    near=[j for j,x in enumerate(fits) if abs(x['score']-selected['score'])<=2*allowance]
    return dict(path=selected['path'],posterior=selected['posterior'],bin_values=values,
        sampling_variance=sampling,total_variance=variance,beta=selected['beta'],tau=tau,centre=centre,centre_mode=centre_mode,
        grid_scores=[dict(beta=x['beta'],score=x['score']) for x in fits],
        near_tied_grid_indices=near if len(near)>1 else [],score_rounding_budget=allowance,
        calibration='Unvalidated conditional RNA boundary model')


def _sum_covering_intervals(values):
    """Sum interval weights at each bin, with exact cancellation at boundaries."""
    n = len(values)
    changes = [0] * (n + 1)
    # Every finite binary64 value is an integer multiple of the least subnormal.
    # Integer range additions preserve the exact sum before rounding each bin.
    unit_exponent = 1074
    denominator = 1 << unit_exponent
    for first, row in enumerate(values):
        for last in range(first, n):
            value = float(row[last])
            if value:
                numerator, divisor = value.as_integer_ratio()
                units = numerator << (unit_exponent - (divisor.bit_length() - 1))
                changes[first] += units
                changes[last + 1] -= units
    total = 0
    sums = np.empty(n)
    for j in range(n):
        total += changes[j]
        sums[j] = total / denominator
    return sums


def _run_probabilities(emission, prior, change):
    """Probability each interval is one complete constant-state run."""
    e, prior = np.asarray(emission, float), np.asarray(prior, float)
    if (e.ndim != 2 or not len(e) or prior.shape != (e.shape[1],)
            or not np.isfinite(e).all() or not np.isfinite(prior).all()
            or np.any(prior <= 0) or not 0 < change < 1):
        raise ValueError('Finite emissions, positive prior and interior change probability required')
    prior = prior/prior.sum()
    n, states = e.shape
    T = np.log(_transition(prior, change))
    forward, backward = np.empty_like(e), np.empty_like(e)
    forward[0] = np.log(prior)+e[0]
    for j in range(1, n):
        forward[j] = e[j]+logsumexp(forward[j-1, :, None]+T, axis=0)
    backward[-1] = 0.
    for j in range(n-2, -1, -1):
        backward[j] = logsumexp(T+e[j+1, None, :]+backward[j+1, None, :], axis=1)
    evidence = float(logsumexp(forward[-1]))
    off = T.copy()
    np.fill_diagonal(off, -np.inf)
    left, right = np.empty_like(e), np.empty_like(e)
    left[0], right[-1] = np.log(prior), 0.
    for j in range(1, n):
        left[j] = logsumexp(forward[j-1, :, None]+off, axis=0)
    for j in range(n-1):
        right[j] = logsumexp(off+e[j+1, None, :]+backward[j+1, None, :], axis=1)
    cumulative = np.vstack([np.zeros(states), np.cumsum(e, axis=0)])
    probabilities = np.zeros((n, n))
    for a in range(n):
        for b in range(a, n):
            logp = left[a]+cumulative[b+1]-cumulative[a]+(b-a)*np.diag(T)+right[b]-evidence
            probabilities[a, b] = float(np.exp(logsumexp(logp)))
    mass = _sum_covering_intervals(probabilities)
    # Bound accumulated log-operation error, its exp propagation and mass sum.
    operations = 128*n*n*states+256
    gamma = operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)
    scale = max(1., abs(evidence), float(np.max(abs(cumulative))), float(np.max(abs(forward))))
    budget = math.expm1(8*gamma*scale)+gamma+n*n*np.nextafter(0., 1.)
    error = float(np.max(abs(mass-1)))
    if not math.isfinite(budget) or error > budget:
        raise ArithmeticError('Unresolved posterior run-mass conservation')
    return probabilities, dict(bins=n, states=states, maximum_mass_error=error,
        mass_rounding_budget=budget, log_evidence=evidence)


class _MassTree:
    """Integer dyadic masses: exact weighted-median ties without float epsilons."""
    def __init__(self, n):
        self.tree = [0]*(n+1)
        self.total = 0

    def add(self, position, value):
        self.total += value
        j = position+1
        while j < len(self.tree):
            self.tree[j] += value
            j += j & -j

    def locate(self, target):
        index, cumulative = 0, 0
        step = 1 << (len(self.tree).bit_length()-1)
        while step:
            nxt = index+step
            if nxt < len(self.tree) and cumulative+self.tree[nxt] < target:
                cumulative += self.tree[nxt]
                index = nxt
            step >>= 1
        return index, cumulative

    def median(self, ordered):
        index, _ = self.locate((self.total+1)//2)
        # The first positive mass at/after index is recovered exactly.
        after = 0
        j = index+1
        while j:
            after += self.tree[j]
            j -= j & -j
        if 2*after == self.total:
            following, _ = self.locate(after+1)
            return float((ordered[index]+ordered[following])/2)
        return float(ordered[index])


def _interval_summaries(values, variance, bins, length, *, weighted_only=False):
    values, variance, bins = np.asarray(values, float), np.asarray(variance, float), np.asarray(bins, int)
    if (values.ndim != 1 or variance.shape != values.shape or bins.shape != values.shape
            or not np.isfinite(values).all() or not np.isfinite(variance).all() or np.any(variance <= 0)
            or not np.array_equal(np.unique(bins), np.arange(length))):
        raise ValueError('Finite gene effects, positive weights and complete local bins required')
    order = np.argsort(values, kind='stable')
    y, v, b = values[order], variance[order], bins[order]
    weights = 1/np.sqrt(v)
    dyadic = [float(w).as_integer_ratio() for w in weights]
    exponent = max(d.bit_length()-1 for _, d in dyadic)
    integers = [p << (exponent-(q.bit_length()-1)) for p, q in dyadic]
    gene_rows = [np.flatnonzero(b == j) for j in range(length)]
    summaries = np.zeros((1 if weighted_only else 3, length, length))
    for a in range(length):
        weighted = _MassTree(len(y))
        uniform = None if weighted_only else _MassTree(len(y))
        selected = []
        for end in range(a, length):
            for j in gene_rows[end]:
                weighted.add(int(j), integers[j])
                if uniform is not None:
                    uniform.add(int(j), 1)
                    selected.append(int(j))
            summaries[0 if weighted_only else 1, a, end] = weighted.median(y)
            if uniform is not None:
                summaries[0, a, end] = uniform.median(y)
                summaries[2, a, end] = (math.fsum(float(y[j]/v[j]) for j in selected)
                                        / math.fsum(float(1/v[j]) for j in selected))
    return summaries


def _average(values, variance, gene_bins, reset, fitted, params, *, weighted_only=False):
    if fitted.get('centre', 0.) != 0.:
        raise ValueError('Only the validated zero-centre working model is supported')
    y, total = np.asarray(fitted['bin_values']), np.asarray(fitted['total_variance'])
    means = fitted['beta']*np.asarray(profile.LOG_DOSE, float)
    emission = -.5*(np.log(2*np.pi*total)[:, None]+(y[:, None]-means)**2/total[:, None])
    starts = np.flatnonzero(reset)
    ends = np.r_[starts[1:], len(y)]
    result = np.empty((1 if weighted_only else 3, len(y)))
    records = []
    for a, b in zip(starts, ends):
        probability, record = _run_probabilities(emission[a:b], profile.STATE_PRIOR, params.p_change)
        keep = (gene_bins >= a)&(gene_bins < b)
        summaries = _interval_summaries(values[keep], variance[keep], gene_bins[keep]-a, b-a,
                                        weighted_only=weighted_only)
        for k, summary in enumerate(summaries):
            result[k, a:b] = _sum_covering_intervals(probability * summary)
        record.update(first_bin=int(a), last_bin=int(b-1))
        records.append(record)
    return result, records


def posterior_rna_profile(gene_log_ratios, gene_expected, gene_bins,
                          reference_log_ratios, reference_expected, blocks, *,
                          reference_id, counts_kind, params=bulk.BULK_PROFILE):
    """Estimate a descriptive bulk RNA profile from prepared gene contrasts.

    Query arrays are genes long; reference arrays are references by genes, in
    the SAME gene order. Expected values are fitted gene exposures, not counts.
    Declare the source as raw or reconstructed; neither acquires calibrated
    count uncertainty from this descriptive function. Apply the same upstream
    preprocessing to query and reference effects. The caller supplies an
    independently chosen reference identity, never selected using DNA scores.

    gene_bins are complete zero-based integers. blocks has one arm/gap label
    per bin; each label must occupy one contiguous interval. At least three
    reference profiles are required for the existing leave-one-out working
    scales. Reference composition and regional transcription remain potential
    explanations for the resulting RNA contrasts.

    The retained normal-centred, reference-weighted Gaussian boundary model
    averages gene medians over all constant-state runs. Its state grid, tau
    floor and variance weights are working smoothing conventions, not DNA
    states, calibrated probabilities, confidence intervals or detection limits.
    For a block of B bins and G genes, interval storage is O(B² + G).
    Run probability and interval accumulation cost O(B²) at fixed state count.
    Weighted medians add O((BG + B²) log(1 + G)) work per block.

    Returns (bin_profile, diagnostics). The optional method reproduced a bulk
    profile regression advantage over pseudobulk; its broader superiority and
    uncertainty gates failed. It changes no existing inference/call defaults.
    """
    if counts_kind not in ('raw', 'reconstructed'):
        raise ValueError('Declare raw or reconstructed source measurements')
    if not isinstance(reference_id, str) or not reference_id.strip():
        raise ValueError('A nonempty explicit reference identity is required')
    values, expected = np.asarray(gene_log_ratios, float), np.asarray(gene_expected, float)
    reference, ref_expected = np.asarray(reference_log_ratios, float), np.asarray(reference_expected, float)
    bins, blocks = np.asarray(gene_bins), np.asarray(blocks)
    if (values.ndim != 1 or not len(values) or expected.shape != values.shape
            or bins.shape != values.shape or bins.dtype.kind not in 'iu'
            or blocks.ndim != 1 or not len(blocks)
            or not np.array_equal(np.unique(bins), np.arange(len(blocks)))
            or reference.ndim != 2 or reference.shape[1] != len(values)
            or len(reference) < 3 or ref_expected.shape != reference.shape
            or not all(np.isfinite(a).all() for a in (values, expected, reference, ref_expected))
            or np.any(expected < 0) or np.any(ref_expected < 0)
            or pd.isna(blocks).any()):
        raise ValueError('Aligned finite gene effects/exposures, three references and complete integer bins/blocks required')
    reset = np.r_[True, blocks[1:] != blocks[:-1]]
    if not pd.Index(blocks[reset]).is_unique:
        raise ValueError('Each arm/gap block label must occupy one contiguous interval')
    if not isinstance(params, profile.ProfileParams):
        raise ValueError('ProfileParams required for the working boundary model')
    if (not np.isfinite(params.tau_floor) or params.tau_floor <= 0
            or not np.isfinite(params.p_change) or not 0 < params.p_change < 1):
        raise ValueError('Positive finite working tau floor and interior change probability required')
    centre = np.median(reference, axis=0)
    spread = (np.median(abs(reference-centre), axis=0)/ndtri(.75))**2
    effect, _, location = _normal_centered_effect(values, reference, bins)
    variance = _gene_variance(spread, expected, 'reference')
    training = _train(reference, ref_expected, bins)
    boundary = _fit(effect, expected, variance, bins, training['reference_variance'], reset, params)
    averaged, records = _average(effect, variance, bins, reset, boundary, params, weighted_only=True)
    diagnostics = dict(method='reference_weighted_posterior_region_median', status='estimated',
        units='natural_log_relative_rna_rate', counts_kind=counts_kind,
        reference_id=reference_id, reference_profiles=len(reference), genes=len(values), bins=len(blocks),
        blocks=int(reset.sum()), centre_rna_log_ratio=location,
        working_parameters=asdict(params), boundary={k: v for k, v in boundary.items()
            if k not in ('path', 'posterior', 'bin_values', 'sampling_variance', 'total_variance')},
        run_mass_diagnostics=records, calibration_status='unvalidated',
        interpretation='Descriptive conditional RNA smoother; no DNA state, clone membership or calibrated interval')
    return averaged[0], diagnostics


def average_rna_profiles(profiles, region_ids, *, reference_ids):
    """Equal mean of independently prepared, identity-aligned RNA profiles.

    Each profile has its own ordered region_ids. All fits must describe exactly
    the same regions and finite support after alignment; missing values cannot
    silently change the estimator. Output follows the first fit's region order.
    reference_ids identify the distinct fitted estimators, not biological n.

    Map EACH fit to the intended common regions before averaging. Averaging
    bins and then taking window medians is a different nonlinear estimator.
    This finite ensemble provides no calibrated uncertainty or DNA states.
    """
    arrays = [np.asarray(a, float) for a in profiles]
    axes = [pd.Index(ids) for ids in region_ids]
    identities = pd.Index(reference_ids)
    if (not arrays or len(arrays) != len(axes) or len(identities) != len(arrays)
            or identities.hasnans or not identities.is_unique
            or any(not str(v).strip() for v in identities)):
        raise ValueError('One unique nonmissing fit identity and region axis per profile required')
    first, aligned = axes[0], []
    for array, axis in zip(arrays, axes):
        if (array.ndim != 1 or not len(array) or len(axis) != len(array)
                or axis.hasnans or not axis.is_unique or len(axis) != len(first)
                or any(not str(v).strip() for v in axis) or np.isinf(array).any()):
            raise ValueError('Nonempty profiles with unique matching region identities required; infinity unsupported')
        index = axis.get_indexer(first)
        if np.any(index < 0):
            raise ValueError('Every fitted profile must cover the identical genomic regions')
        aligned.append(array[index])
    aligned = np.asarray(aligned)
    support = np.isfinite(aligned)
    if not np.all(support == support[:1]):
        raise ValueError('All fitted profiles must have identical finite support')
    result = np.full(len(first), np.nan)
    for j in np.flatnonzero(support[0]):
        result[j] = math.fsum(map(float, aligned[:, j]))/len(arrays)
    return result, dict(method='equal_reference_profile_mean', reference_ids=identities.tolist(),
        region_ids=first.tolist(), fitted_estimators=len(arrays), finite_regions=int(support[0].sum()),
        units='natural_log_relative_rna_rate', calibration_status='unvalidated',
        interpretation='Finite reference ensemble; fit count is not biological replication or a posterior sample size')
