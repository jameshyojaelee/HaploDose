"""Observe or explicitly apply research gene-pool alignment; no public default."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys
import weakref

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import binom
from haplodose.clones import SpecimenModel
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def decision(ya, yb, gbin, mask, seg, min_genes=8, alpha=.01):
    if ya.sum() <= 0 or yb.sum() <= 0:
        return dict(available=False, reason='zero_pool')
    r = np.log((ya+.5)/ya.sum())-np.log((yb+.5)/yb.sum())
    informative = ya+yb >= 10
    if informative.sum() < min_genes:
        return dict(available=True, passed=True, reason='insufficient_informative_genes')
    r -= np.median(r[informative])
    selected = informative & (gbin >= seg['first_bin']) & (gbin <= seg['last_bin']) & mask[np.clip(gbin, 0, len(mask)-1)]
    n = int(selected.sum())
    k = int(sum(np.sign(r[selected]) == (1 if seg['direction'] == 'higher' else -1)))
    return dict(available=True, passed=bool(n < min_genes or binom.sf(k-1, n, .5) < alpha), n=n, k=k,
        p=None if n < min_genes else float(binom.sf(k-1, n, .5)),
        reason='insufficient_segment_genes' if n < min_genes else 'sign_test')


def effective_counts(ya, yb, weights):
    """Use the regional conditional-Poisson weight unit, not a new count law.

    c=sum(w²)/sum(w); S/c and L/c preserve weighted population rates, restore
    the raw-count rule for binary weights, and remove common weight units.
    This does not calibrate gene signs or account for NB/selection uncertainty.
    """
    W = np.asarray(weights, np.float64)
    if W.ndim != 2 or W.shape[1] != 2 or not np.isfinite(W).all() or np.any(W < 0):
        raise ValueError('Expected finite nonnegative two-pool weights')
    mass = W.sum(axis=0)
    square = (W*W).sum(axis=0)
    if np.any(mass <= 0) or np.any(square <= 0) or not np.isfinite(square).all():
        raise ValueError('Positive finite pool weight units required')
    c = square/mass
    return np.asarray(ya, np.float64)/c[0], np.asarray(yb, np.float64)/c[1], c


@contextmanager
def observe_weights(trace, directory, *, apply=False, weight_units='raw'):
    if weight_units not in ('raw', 'effective'):
        raise ValueError('Unknown research weight units')
    M = SpecimenModel
    pooled, veto = M.pooled, M.gene_consistent
    cache, models = weakref.WeakKeyDictionary(), weakref.WeakKeyDictionary()

    def pools(self, *args, **kwargs):
        caller = sys._getframe(1)
        if caller.f_code.co_name == 'matched_differences' and caller.f_locals.get('self') is self:
            loc = caller.f_locals
            if 'W' not in loc or args[0] is not loc['W']:
                raise RuntimeError('Unexpected matched pooling identity')
            registry = cache.setdefault(self, {})
            for i, (a, b) in enumerate(loc['pairs']):
                key = (np.asarray(loc['labels'], np.int64).tobytes(), int(a), int(b))
                registry[key] = loc['W'][:, 2*i:2*i+2].toarray()
        return pooled(self, *args, **kwargs)

    def gene(self, seg, labels, a, b, min_genes=8, alpha=.01):
        result = veto(self, seg, labels, a, b, min_genes=min_genes, alpha=alpha)
        if self.genes is None:
            return result
        sub, X, gbin = self.genes
        if self not in models:
            index = len(list(directory.glob('model_*.npz')))
            path = directory/f'model_{index}.npz'
            X = sparse.csr_matrix(X)
            seeds = [int(part.removeprefix('seed_')) for part in Path(self.rows.Y.filename).parts if part.startswith('seed_')]
            assert len(seeds) == 1
            np.savez_compressed(path, data=X.data, indices=X.indices, indptr=X.indptr, shape=np.array(X.shape),
                sub=sub, gbin=gbin, rows=self.rows.rows, mask=self.test_mask, seed=np.array(seeds[0]))
            models[self] = str(path)
        pairlabels = np.asarray(labels)[sub]
        ya = np.asarray(X[pairlabels == a].sum(axis=0)).ravel()
        yb = np.asarray(X[pairlabels == b].sum(axis=0)).ravel()
        key = (np.asarray(labels, np.int64).tobytes(), int(a), int(b))
        W = cache.get(self, {}).get(key)
        hard = decision(ya, yb, gbin, self.test_mask, seg, min_genes, alpha)
        # Original early cell-support abstention occurs before gene counts.
        early = min(sum(pairlabels == a), sum(pairlabels == b)) < 10
        assert early or hard['passed'] == bool(result)
        weighted = dict(available=False, reason='no_corresponding_matched_pool')
        saved = dict(labels=labels, hard_a=ya, hard_b=yb)
        if W is not None:
            represented = np.zeros(self.rows.n, bool); represented[sub] = True
            assert np.all(W >= 0) and np.isfinite(W).all()
            if np.any(W[~represented] != 0):
                weighted = dict(available=False, reason='gene_subset_does_not_span_regional_pool')
            else:
                wa, wb = W[sub, 0], W[sub, 1]
                xa, xb = np.asarray(X.T@wa).ravel(), np.asarray(X.T@wb).ravel()
                saved.update(weights=W, weighted_a=xa, weighted_b=xb)
                if weight_units == 'effective':
                    xa, xb, c = effective_counts(xa, xb, W[sub])
                    saved.update(effective_a=xa, effective_b=xb, weight_scale=c)
                weighted = decision(xa, xb, gbin, self.test_mask, seg, min_genes, alpha)
        path = directory/f'veto_{len(trace):04d}.npz'
        np.savez_compressed(path, **saved)
        chosen = bool(weighted['passed']) if apply and not early and weighted['available'] else bool(result)
        trace.append(dict(index=len(trace), model=models[self], arrays=str(path), seg=dict(seg),
            a=int(a), b=int(b), min_genes=min_genes, alpha=alpha, hard=hard, weighted=weighted,
            early_cell_abstention=bool(early), original_passed=bool(result), applied_passed=chosen,
            weight_units=weight_units,
            caller=sys._getframe(1).f_back.f_code.co_name))
        return chosen

    M.pooled, M.gene_consistent = pools, gene
    try:
        yield
    finally:
        M.pooled, M.gene_consistent = pooled, veto


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--apply', action='store_true', help='Research-only aligned-pool veto, no calibration claim')
    p.add_argument('--weight-units', choices=('raw','effective'), default='raw')
    a = p.parse_args(); a.out.mkdir(exist_ok=False)
    directory = a.out/'arrays'; directory.mkdir()
    root = Path(__file__).parent
    spec = importlib.util.spec_from_file_location('fixed_gene_runner', root/'gene_thinning_recovery.py')
    runner = importlib.util.module_from_spec(spec); spec.loader.exec_module(runner)
    rec = start_run(a.out, dict(script=Path(__file__), runner=root/'gene_thinning_recovery.py',
        hook=root/'gene_thinning_hook.py', stage=root/'membership_weighted_stages.py',
        design=a.task/'DESIGN.md'), dict(apply=a.apply, weight_units=a.weight_units), dict(scope='Fixed-fold population diagnostic or explicit research intervention'), mode='gene_veto_weights')
    trace, argv = [], sys.argv
    try:
        sys.argv = [str(runner.__file__), '--out', str(a.out/'replay'), '--task',
            str(a.task.parent/'gene_thinning_20261109'), '--mode', 'gene', '--stage', 'both']
        if a.apply:
            assert json.loads((a.task.parent/'gene_veto_weights_20261109/audit_v1/report.json').read_text())['all_checks_passed']
        with observe_weights(trace, directory, apply=a.apply, weight_units=a.weight_units):
            runner.main()
        previous = a.task.parent/'gene_thinning_weighted_20261109/both_v1'
        current = a.out/'replay'
        old = json.loads((previous/'panel/summary.json').read_text())
        new = json.loads((current/'panel/summary.json').read_text())
        if not a.apply:
            for before, after in zip(old['cases'], new['cases']):
                assert before['seed'] == after['seed'] and before['criteria'] == after['criteria']
                assert before['original_assertion_passed'] == after['original_assertion_passed']
                for name in ('cells.tsv.gz', 'differences.tsv'):
                    pd.testing.assert_frame_equal(pd.read_csv(previous/f'panel/seed_{before["seed"]}/fit'/name, sep='\t'),
                        pd.read_csv(current/f'panel/seed_{before["seed"]}/fit'/name, sep='\t'), check_exact=True)
        write_json(dict(apply=a.apply, all_fit_tables_and30_criteria_identical=not a.apply, calls=trace, biological_specimens=0,
            production_changed=False, scientific_criteria_changed=False), a.out/'report.json')
    except BaseException as error:
        write_json(trace, a.out/'partial_trace.json')
        finish_run(a.out, rec, error=error)
        raise
    finally:
        sys.argv = argv
    finish_run(a.out, rec, backend='CPU fixed known-noise workflow; veto returns unchanged')


if __name__ == '__main__':
    main()
