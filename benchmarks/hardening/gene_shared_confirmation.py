"""Research exact cell-shared membership on gene or moment-bin NB evidence."""
from contextlib import contextmanager
import hashlib
import importlib.util
from pathlib import Path
import sys

import numpy as np
from haplodose import clones, profile
from haplodose.io import write_json


def load(name):
    spec = importlib.util.spec_from_file_location('gene_shared_'+name, Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gene_model(O, E, r, q, gene_bins, effects, prior, change, *, left_neutral,
               right_neutral, mixture_proportion):
    """Sum gene log likelihoods before introducing one chain state per bin."""
    engine = load('shared_membership_hmm')
    bins = np.asarray(gene_bins, int)
    length = len(mixture_proportion)
    if bins.shape != (np.shape(O)[1],) or np.any(bins < 0) or np.any(bins >= length):
        raise ValueError('Invalid gene-to-region mapping')
    if set(bins) != set(range(length)):
        raise ValueError('Every region bin requires identified gene composition')
    model = engine.build(O, E, r, q, effects, prior, change,
        left_neutral=left_neutral, right_neutral=right_neutral,
        mixture_proportion=np.asarray(mixture_proportion)[bins])
    ll = np.zeros((2, len(q), length, 3))
    for j in range(length):
        ll[:, :, j, :] = model['ll'][:, :, bins == j, :].sum(axis=2)
    model['ll'] = ll
    model['prior_normalizer'] = engine.chain(np.zeros((length, 3)), model['logT'],
        model['initial'], model['terminal'])[1]
    return model


def thin_matching(counts, size, weights, seed):
    """NB count thinning, retaining every positive-weight cell's identity."""
    counts = np.asarray(counts)
    size = np.broadcast_to(np.asarray(size, float), counts.shape)
    weights = np.asarray(weights, float)
    if (counts.ndim != 2 or weights.shape != (len(counts),) or not np.isfinite(counts).all()
            or np.any(counts < 0) or np.any(counts != np.floor(counts))
            or not np.isfinite(size).all() or np.any(size <= 0)
            or not np.isfinite(weights).all() or np.any(weights < 0) or np.any(weights > 1)):
        raise ValueError('Integer counts, positive NB sizes and matching weights in[0,1] required')
    result = np.zeros(counts.shape, dtype=np.int64)
    one, interior = weights == 1, (weights > 0)&(weights < 1)
    result[one] = counts[one].astype(np.int64)
    rng = np.random.default_rng(seed)
    w = weights[interior, None]
    probability = rng.beta(size[interior]*w, size[interior]*(1-w))
    result[interior] = rng.binomial(counts[interior].astype(np.int64), probability)
    if np.any(result > counts):
        raise ArithmeticError('Count matching violates conservation')
    return result


@contextmanager
def confirmation(registry, mode, trace, directory):
    if mode not in ('control', 'moment_bin', 'gene', 'count_bin', 'count_gene', 'depth_conditioned', 'depth_gene_score'):
        raise ValueError(mode)
    use_conditioning=mode in ('depth_conditioned','depth_gene_score')
    M = clones.SpecimenModel
    original_fit, event_code = clones._fit_event_mixture, M.event_mixtures.__code__
    captured, composition_cache = {}, {}
    arrays = load('gene_event_proposal')
    engine, matching = load('shared_membership_hmm'), load('shared_membership_stages')
    active = None
    native_matched = M.matched_differences

    def conditioned_matched(self, labels, pairs, *args, **kwargs):
        if active is None or active[0] is not self:
            return native_matched(self, labels, pairs, *args, **kwargs)
        if args or len(pairs)!=1:
            raise ValueError('One binary depth-conditioned contrast required')
        # The new likelihood uses unit-weight common-support cells. Its support
        # must describe those actual observations, not the superseded thinning
        # margins. The minimum itself and all downstream vetoes stay unchanged.
        options=dict(kwargs);minimum=options.pop('min_effective_cells',self.params.min_cells)
        if minimum!=self.params.min_cells:raise ValueError('Unchanged minimum support required')
        result=native_matched(self,labels,pairs,min_effective_cells=None,**options)
        fitted=active[1]['conditioned_parameters'];groups=fitted['contributing_groups']
        for row in result:
            row.update(effective_cells=groups,depth_adjustment='discovery_conditioned_raw_counts')
            if not fitted['identified'] or min(groups)<minimum:
                row.update(segments=[],testable=False,untestable_reason='insufficient_conditioned_common_support')
        return result

    def fit(O, E, r, step, *, maxiter=1000):
        result = original_fit(O, E, r, step, maxiter=maxiter)
        frame = sys._getframe(1)
        if frame.f_code is event_code:
            loc = frame.f_locals
            model, cells = loc['self'], np.asarray(loc['cells'], int)
            cols = loc['cols'][loc['part']]
            source = f"event_chrom{loc['c']}_bins{cols[0]}-{cols[-1]}"
            captured[(id(model.rows), cells.tobytes(), source)] = dict(
                E=E.copy(), q=result['responsibility'].copy(), cols=cols.copy())
        return result

    clones._fit_event_mixture = fit
    if use_conditioning:M.matched_differences=conditioned_matched
    try:
        # Load after installing fit, so its existing caller-identity observer
        # observes this exact wrapper and the native event frame consistently.
        with load('membership_weighted_stages').weighted_stages('both', trace):
            score, difference = M._score_candidate, M.differences
            original_gene=M.gene_consistent

            def gene_concordance(self, segment, labels, a, b, min_genes=8, alpha=.01):
                old=original_gene(self,segment,labels,a,b,min_genes=min_genes,alpha=alpha)
                if active is None or active[0] is not self:return old
                if (a,b)!=(0,1) or min_genes!=8 or alpha!=.01:
                    raise ValueError('Original binary gene-support criteria required')
                _,item,cells,g,entry=active
                np.testing.assert_array_equal(np.asarray(labels)[cells],g)
                if self.genes is None:raise RuntimeError('Gene-score mode requires gene measurements')
                observation=item['conditioned_gene_observation'];parameters=item['conditioned_parameters']
                choose=np.asarray(parameters['selected'],int);counts=observation['counts'][choose]
                sizes=np.broadcast_to(observation['size'],observation['counts'].shape)[choose]
                expected=observation['expected'][choose]*np.exp(np.asarray(parameters['offsets'])[choose,None])
                result=load('conditioned_gene_score').concordance(counts,expected,sizes,item['q'][choose],
                    np.asarray(parameters['strata'],int)[choose],observation['bins'],
                    self.test_mask[observation['bins']],segment,min_genes=min_genes,alpha=alpha)
                records=entry.setdefault('gene_score_checks',[])
                path=directory/f'gene_score_{len(trace):05d}_{len(records):02d}.json'
                if path.exists():raise RuntimeError('Gene-score record already exists')
                write_json(dict(segment=segment,proposal_lock=item['conditioned_lock'],
                    original_applied_veto_passed=bool(old),new_gene_score=result),path)
                records.append(dict(path=str(path),segment=dict(segment),original_applied_veto_passed=bool(old),
                    new_passed=result['passed'],n=result['n'],k=result['k'],nominal_p=result['p']))
                return bool(result['passed'])

            def scored(self, k, cells, g, source):
                nonlocal active
                if mode == 'control' or not source.startswith('event_'):
                    return score(self, k, cells, g, source)
                cells = np.asarray(cells, int)
                item = captured[(id(self.rows), cells.tobytes(), source)]
                if not np.array_equal(item['q'] > .5, g.astype(bool)):
                    raise RuntimeError('Membership/order mismatch')
                entry = dict(stage='gene_shared_confirmation', mode=mode,
                    source=source, case_seed=registry[id(self.rows)]['case_seed'], calls=[])
                trace.append(entry)
                if 3**len(item['cols']) > 1000000:
                    entry['fallback'] = 'declared_exact_state_cap'
                    result = score(self, k, cells, g, source)
                else:
                    if use_conditioning:
                        info=registry[id(self.rows)];cols=item['cols'];X=info['prop'][cells].toarray()
                        fraction,size,support=arrays.composition(X,info['gbin'],
                            info['resolved_gene_size']*(1-info['share']),self.J)
                        composition_cache[(id(self.rows),cells.tobytes())]=(fraction,size,support)
                        if not support[cols].all():
                            entry.update(accepted=False,reason='unidentified_proposal_gene_composition');return None
                        O,Eg,rg,genes,local_bin=arrays.region_arrays(X,info['gbin'],fraction,
                            info['resolved_gene_size']*info['share'],cols,item['E'])
                        outside=~np.isin(info['gbin'],cols)
                        if not outside.any():
                            entry.update(accepted=False,reason='no_outside_event_proposal_depth');return None
                        depth=np.log(np.maximum(X[:,outside].sum(axis=1),1.))
                        parameters=load('depth_conditioned_genes').prepare(O,Eg,rg,item['q'],depth,min_cells=self.params.min_cells)
                        lock=directory/f'proposal_lock_{len(trace):05d}.json'
                        if lock.exists():raise RuntimeError('Proposal lock already exists')
                        write_json(dict(parameters=parameters,genes=genes.tolist(),
                            rows=self.rows.rows[cells].tolist(),source=source,
                            confirmation_pooling_started=False,confirmation_count_arrays_opened=False),lock)
                        item.update(conditioned_parameters=parameters,conditioned_genes=genes,conditioned_lock=str(lock))
                        entry.update(proposal_lock=str(lock),contributing_groups=parameters['contributing_groups'])
                    old = active
                    active = (self, item, cells, np.asarray(g), entry)
                    try:
                        result = score(self, k, cells, g, source)
                    finally:
                        active = old
                entry['accepted'] = result is not None
                return result

            def differences(self, *args, **kwargs):
                if active is None or active[0] is not self:
                    return difference(self, *args, **kwargs)
                _, item, cells, g, entry = active
                caller = sys._getframe(1)
                matched = caller.f_code.co_name == 'matched_differences' and 'W' in caller.f_locals
                context = 'depth_matched' if matched else 'initial'
                saved = {}
                def observe(frame, event, value):
                    if event == 'return' and frame.f_code is profile.difference_segments.__code__:
                        saved.update({k:frame.f_locals[k].copy() for k in ('d', 'var', 'idx')})
                previous = sys.getprofile()
                sys.setprofile(observe)
                try:
                    old = difference(self, *args, **kwargs)
                finally:
                    sys.setprofile(previous)
                if len(old) != 1 or saved['d'].shape[0] != 1:
                    raise RuntimeError('Binary comparison required')
                info, cols, q = registry[id(self.rows)], item['cols'], item['q']
                if np.any(np.diff(cols) != 1) or np.any(self.reset[cols[1:]]):
                    raise ValueError('Event crosses genomic discontinuity')
                key = (id(self.rows), cells.tobytes())
                if key not in composition_cache:
                    X = info['prop'][cells].toarray()
                    fraction, size, support = arrays.composition(X, info['gbin'],
                        info['resolved_gene_size']*(1-info['share']), self.J)
                    composition_cache[key] = (fraction, size, support)
                fraction, size, support = composition_cache[key]
                call = dict(context=context, original_segments=old[0]['segments'])
                entry['calls'].append(call)
                def abstain(reason):
                    call.update(reason=reason, exact_segments=[])
                    return [dict(old[0], segments=[], testable=False, untestable_reason=reason)]
                if not self.test_mask[cols].all():
                    return abstain('uncallable_bin_in_event')
                if not support[cols].all():
                    return abstain('unidentified_gene_composition')
                choose, sizes, matching_weights = np.arange(len(cells)), None, None
                if matched and mode.startswith('count_'):
                    W = caller.f_locals['W']
                    if W.shape != (self.rows.n, 2):
                        raise RuntimeError('One matched contrast required')
                    matching_weights = np.asarray(W[cells].sum(axis=1)).ravel()
                    if np.any(matching_weights > 1):
                        raise ValueError('Executed matching weight exceeds one')
                    choose = np.flatnonzero(matching_weights > 0)
                elif matched and not use_conditioning:
                    choose, sizes = matching.balanced_subset(self.x[cells], g, self.rows.rows[cells])
                call.update(selected_cells=choose.tolist(), per_stratum_per_group=sizes)
                if min(np.bincount(g[choose], minlength=2)) < self.params.min_cells:
                    return abstain('insufficient_depth_overlap')
                scale = self.fold_scale['test']/self.fold_scale['prop']
                E = item['E']*scale
                pi = (E.T@q)/E.sum(axis=0)
                effects = profile.DIFF_LEVELS*self.beta*np.log(1.5)
                options = dict(left_neutral=bool(not self.reset[cols[0]]),
                    right_neutral=bool(cols[-1]+1 < self.J and not self.reset[cols[-1]+1]),
                    mixture_proportion=pi)
                if mode in ('gene', 'count_gene', 'depth_conditioned', 'depth_gene_score'):
                    X = (info['X'][cells]-info['prop'][cells]).toarray()
                    O, Eg, rg, genes, local_bin = arrays.region_arrays(X, info['gbin'],
                        fraction, info['resolved_gene_size']*(1-info['share']), cols, E)
                    if np.any((fraction == 0) & (X.sum(axis=0) > 0) & np.isin(info['gbin'], cols)):
                        return abstain('unseen_gene_has_confirmation_counts')
                    summed_O = np.stack([O[:, local_bin == j].sum(axis=1) for j in range(len(cols))], axis=1)
                    summed_E = np.stack([Eg[:, local_bin == j].sum(axis=1) for j in range(len(cols))], axis=1)
                    np.testing.assert_array_equal(summed_O, self.rows.take(cells, self.test_fold)[:, cols])
                    operations = 4*len(genes)+32
                    error_bound = operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)*max(1., float(np.max(abs(E))))
                    if np.max(abs(summed_E-E)) > error_bound:
                        raise ArithmeticError('Gene expectations do not conserve bin expectations')
                else:
                    O = self.rows.take(cells, self.test_fold)[:, cols]
                    Eg, rg, genes, local_bin = E, size[cols], np.array([], int), np.arange(len(cols))
                thinning_seed = None
                if matching_weights is not None:
                    tag = f"HaploDose-count-matching-v1:{info['case_seed']}:{entry['source']}:{self.params.seed}".encode()
                    tag += np.asarray(self.rows.rows[cells], dtype='<i8').tobytes()
                    thinning_seed = int.from_bytes(hashlib.sha256(tag).digest()[:8], 'little')
                    O = thin_matching(O, rg, matching_weights, thinning_seed)
                    Eg = Eg*matching_weights[:, None]
                    rg = np.broadcast_to(rg, O.shape)*matching_weights[:, None]
                fitted_size = rg[choose] if np.ndim(rg) == 2 else rg
                if use_conditioning:
                    np.testing.assert_array_equal(genes,item['conditioned_genes'])
                    conditioned=item['conditioned_parameters']
                    if not conditioned['identified']:return abstain('insufficient_conditioned_common_support')
                    choose=np.asarray(conditioned['selected'],int)
                    logT=np.log(engine.transition(profile.DIFF_PRIOR,self.diff_pp.p_change))
                    initial=logT[1].copy() if options['left_neutral'] else np.log(profile.DIFF_PRIOR)
                    terminal=logT[:,1].copy() if options['right_neutral'] else np.zeros(3)
                    model=load('depth_conditioned_genes').model(O,Eg,rg,q,local_bin,effects,profile.DIFF_PRIOR,
                        self.diff_pp.p_change,conditioned,initial,terminal)
                    item['conditioned_gene_observation']=dict(counts=O,expected=Eg,size=rg,bins=cols[local_bin])
                    call.update(proposal_lock=item['conditioned_lock'],observation_model='depth_conditioned_raw_gene_counts',
                        contributing_groups=conditioned['contributing_groups'])
                elif mode in ('gene', 'count_gene'):
                    model = gene_model(O[choose], Eg[choose], fitted_size, q[choose], local_bin,
                        effects, profile.DIFF_PRIOR, self.diff_pp.p_change, **options)
                else:
                    model = engine.build(O[choose], Eg[choose], fitted_size, q[choose], effects,
                        profile.DIFF_PRIOR, self.diff_pp.p_change, **options)
                fitted = engine.enumerate_states(model)
                reset = np.r_[True, np.zeros(len(cols)-1, bool)]
                # Callable-mask exclusions terminate cores rather than bridging
                # missing evidence. They do not receive artificial posteriors.
                valid = self.test_mask[cols]
                reset[1:] |= ~valid[:-1] | ~valid[1:]
                posterior, path = fitted['posterior'].copy(), fitted['path'].copy()
                posterior[~valid] = [0., 1., 0.]
                path[~valid] = 1
                segments = []
                for a, b, s in profile.credible_difference_runs(path, posterior, reset, self.diff_pp):
                    first, last = int(cols[a]), int(cols[b])
                    use = (saved['idx'] >= first) & (saved['idx'] <= last)
                    w = 1/saved['var'][0, use]
                    mean = float(np.sum(w*saved['d'][0, use])/w.sum())
                    prob = posterior[a:b+1, s]
                    segments.append(dict(first_bin=first, last_bin=last, n_test_bins=b-a+1,
                        direction='lower' if s == 2 else 'higher', mean_diff=mean,
                        z=float(mean*np.sqrt(w.sum())), posterior=float(prob.mean()),
                        min_bin_posterior=float(prob.min()), posterior_rule='every_retained_bin'))
                payload = directory/f'comparison_{sum(len(t.get("calls", [])) for t in trace):05d}.npz'
                if payload.exists():
                    raise RuntimeError('Comparison output already exists')
                np.savez_compressed(payload, counts=O, expected=Eg, size=rg, cells=cells,
                    global_rows=self.rows.rows[cells], selected=choose, q=q, cols=cols,
                    genes=genes, local_bin=local_bin, fractions=fraction, mixture_proportion=pi,
                    effects=effects, prior=profile.DIFF_PRIOR, change=self.diff_pp.p_change,
                    matching_weights=np.ones(len(cells)) if matching_weights is None else matching_weights,
                    ll=model['ll'], initial=model['initial'], terminal=model['terminal'],
                    posterior=fitted['posterior'], path=fitted['path'], membership=fitted['membership'])
                call.update(reason='exact_conditional', exact_segments=segments,
                    arrays=str(payload), states_enumerated=fitted['states_enumerated'],
                    thinning_seed=thinning_seed, represented_cells=len(choose),
                    log_evidence=fitted['log_evidence'], max_nonzero_posterior=float(fitted['posterior'][:, [0,2]].max()))
                return [dict(old[0], segments=segments)]

            M._score_candidate, M.differences = scored, differences
            if mode=='depth_gene_score':M.gene_consistent=gene_concordance
            try:
                yield
            finally:
                M._score_candidate, M.differences = score, difference
                if mode=='depth_gene_score':M.gene_consistent=original_gene
    finally:
        clones._fit_event_mixture = original_fit
        if use_conditioning:M.matched_differences=native_matched
