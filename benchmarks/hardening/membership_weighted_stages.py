"""Experimental stage ablations, never installed or enabled by default.

Weights are classification margins, not latent-component deconvolution. Fitted
weights and nuisance parameters remain selected; no independence claim is made.
The unchanged complete workflow must establish error and power empirically.
"""
from contextlib import contextmanager
import sys
import numpy as np
from scipy import sparse
from haplodose.clones import SpecimenModel, _fit_event_mixture


@contextmanager
def weighted_stages(stage, trace, *, record_memberships=False, candidate_filter=None):
    if stage not in {'hard', 'confirmation', 'merge', 'both'}:
        raise ValueError(stage)
    if sys.getprofile() is not None:
        raise RuntimeError('Cannot nest event probability profiling')
    original_event = SpecimenModel.event_mixtures
    original_score = SpecimenModel._score_candidate
    original_matched = SpecimenModel.matched_differences
    cache = {}

    def events(self, *args, **kwargs):
        if stage not in {'confirmation', 'both'}:
            return original_event(self, *args, **kwargs)
        if self.alleles is not None or self.graph is not None:
            raise ValueError('Experimental weights support expression-only nonspatial fits')
        cells = kwargs.get('cells', args[5] if len(args) > 5 else None)
        captured = {}
        def observe(frame, event, value):
            if frame.f_code is _fit_event_mixture.__code__ and event == 'return' and isinstance(value, dict):
                parent = frame.f_back
                if parent.f_code is original_event.__code__:
                    loc = parent.f_locals
                    part, cols = loc['part'], loc['cols']
                    src = f"event_chrom{loc['c']}_bins{cols[part[0]]}-{cols[part[-1]]}"
                    captured[src] = dict(q=value['responsibility'].copy())
                    if candidate_filter is not None:
                        captured[src].update(O=loc['O'].copy(), E=loc['E'].copy(), r=loc['r_nb'].copy())
        previous_profile = sys.getprofile()
        sys.setprofile(observe)
        try:
            result = original_event(self, *args, **kwargs)
        finally:
            sys.setprofile(previous_profile)
        if cells is None:
            raise RuntimeError('Explicit cell identities required')
        # Retain only returned candidates; each parent proposal call replaces its cache.
        key = (id(self), np.asarray(cells).tobytes())
        cache[key] = {}
        for g, src in result:
            q = captured[src]['q']
            if q.shape != g.shape or not np.array_equal(q > .5, g.astype(bool)):
                raise RuntimeError('Captured probability identity/order mismatch')
            cache[key][src] = captured[src]
        return result

    def score(self, k, cells, g, source):
        item = cache.get((id(self), np.asarray(cells).tobytes()), {}).get(source)
        q = None if item is None else item['q']
        if stage not in {'confirmation', 'both'} or q is None:
            if stage in {'confirmation', 'both'} and source.startswith('event_'):
                raise RuntimeError('Event probabilities missing; refusing silent hard fallback')
            return original_score(self, k, cells, g, source)
        if self.Ak is not None or g.max() != 1:
            raise ValueError('Only binary expression proposals supported')
        w = np.abs(2*q-1)
        W = sparse.csr_matrix((w, (cells, g)), shape=(self.rows.n, 2))
        Y, L, XL = self.pooled(W, self.test_fold)
        Q = self.rows.weighted_sum(W.multiply(W), self.test_fold)
        sw = np.asarray(W.sum(axis=0)).ravel()
        sw2 = np.asarray(W.multiply(W).sum(axis=0)).ravel()
        c = np.divide(sw2, sw, out=np.zeros(2), where=sw>0)
        first = self.differences(Y[:1], L[:1], XL[:1], Y[1:], L[1:], XL[1:],
                                count_variance=(Q[:1], Q[1:]), weight_scale=(c[:1], c[1:]))[0]
        entry = dict(stage='confirmation', source=source, parent=int(k), cells=len(cells),
                     specimen_rows=self.rows.rows.tolist(), proposal_cells=np.asarray(cells).tolist(),
                     mean_margin=float(w.mean()), initial_segments=len(first['segments']),
                     passed=False)
        trace.append(entry)
        if not first['segments']:
            entry['reason'] = 'initial_weighted_hmm'
            return None
        labels = np.full(self.rows.n, -1); labels[cells] = g
        weights = np.zeros(self.rows.n); weights[cells] = w
        result = original_matched(self, labels, [(0,1)], cell_weights=weights,
                                  min_effective_cells=self.params.min_cells)[0]
        supported = self.supported(result['segments'], None,None,None,None, labels=labels, a=0,b=1)
        value = sum(abs(s['z']) for s in supported)
        entry.update(effective_cells=result['effective_cells'], matched_segments=len(result['segments']),
                     supported_segments=len(supported), score=float(value), passed=bool(value>0),
                     reason='accepted' if value>0 else 'weighted_match_or_gene_veto')
        if value<=0:
            return None
        if candidate_filter is not None and not candidate_filter(self, cells, g, source, item):
            entry.update(passed=False, reason='continuous_state_alternative')
            return None
        # Preserve the original full-count hard-pool initialization. This ablates
        # confirmation only and does not smuggle a profile-estimation change in.
        hard = sparse.csr_matrix((np.ones(len(cells)), (cells,g)), shape=(self.rows.n,2))
        Ym,Lm,XLm = self.pooled(hard)
        return value, [[0],[1]], Ym,Lm,XLm,None

    def matched(self, labels, pairs, *args, **kwargs):
        caller = sys._getframe(1)
        if record_memberships and caller.f_code is SpecimenModel.merge_unsupported.__code__:
            source = getattr(self, '_research_merge_assignment_source', 'proposal' if self.thin else 'full')
            entry = dict(stage='merge_membership', pairs=[list(pair) for pair in pairs],
                specimen_rows=self.rows.rows.tolist(), labels=np.asarray(labels).tolist(),
                full_responsibility=caller.f_locals['R'].tolist(),
                doses=np.asarray(caller.f_locals['doses']).tolist(),
                mixture_weights=np.asarray(caller.f_locals['pi']).tolist(), assignment_source=source)
            key = 'proposal_responsibility' if source == 'proposal' else 'merge_responsibility'
            entry[key] = caller.f_locals['assignment_R'].tolist()
            trace.append(entry)
        if stage not in {'merge', 'both'} or caller.f_code is not SpecimenModel.merge_unsupported.__code__:
            return original_matched(self, labels, pairs, *args, **kwargs)
        if args or kwargs or self.alleles is not None or self.graph is not None:
            raise ValueError('Unsupported weighted-merge configuration')
        R = caller.f_locals['assignment_R']
        results=[]
        for a,b in pairs:
            w = np.abs(R[:,a]-R[:,b])
            r = original_matched(self,labels,[(a,b)],cell_weights=w,
                                 min_effective_cells=self.params.min_cells)[0]
            trace.append(dict(stage='merge', a=int(a),b=int(b),
                              specimen_rows=self.rows.rows.tolist(),
                              mean_margin=float(w[(labels==a)|(labels==b)].mean()),
                              effective_cells=r['effective_cells'], segments=len(r['segments']),
                              testable=bool(r.get('testable',True))))
            results.append(r)
        return results

    SpecimenModel.event_mixtures=events
    SpecimenModel._score_candidate=score
    SpecimenModel.matched_differences=matched
    try:
        yield
    finally:
        SpecimenModel.event_mixtures=original_event
        SpecimenModel._score_candidate=original_score
        SpecimenModel.matched_differences=original_matched
        cache.clear()
