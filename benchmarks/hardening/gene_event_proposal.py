"""Research gene-versus-bin event likelihoods on the same discovery counts."""
from contextlib import contextmanager
import sys

import numpy as np
from haplodose import clones


def composition(counts, gene_bins, gene_sizes, n_bins):
    """Within-bin gene proportions fitted to discovery counts of one parent.

    All-zero genes have zero fitted rate and provide no fitted likelihood term.
    A completely unobserved bin has unresolved gene composition. No positive
    floor, oracle composition or confirmation count is supplied to fill it.
    """
    counts, gene_bins = np.asarray(counts, float), np.asarray(gene_bins, int)
    sizes = np.broadcast_to(np.asarray(gene_sizes, float), (len(gene_bins),))
    if (counts.ndim != 2 or counts.shape[1] != len(gene_bins) or not np.isfinite(counts).all()
            or np.any(counts < 0) or np.any(gene_bins < 0) or np.any(gene_bins >= n_bins)
            or not np.isfinite(sizes).all() or np.any(sizes <= 0)):
        raise ValueError('Aligned nonnegative counts, valid bins and positive gene sizes required')
    totals = counts.sum(axis=0)
    bin_totals = np.bincount(gene_bins, weights=totals, minlength=n_bins)
    fractions = np.divide(totals, bin_totals[gene_bins], out=np.zeros(len(gene_bins)),
                          where=bin_totals[gene_bins] > 0)
    inverse_size = np.bincount(gene_bins, weights=fractions**2/sizes, minlength=n_bins)
    bin_sizes = np.divide(1., inverse_size, out=np.full(n_bins, np.nan), where=inverse_size > 0)
    return fractions, bin_sizes, bin_totals > 0


def region_arrays(counts, gene_bins, fractions, gene_sizes, selected_bins, expected_bins):
    selected_bins, gene_bins = np.asarray(selected_bins, int), np.asarray(gene_bins, int)
    if len(np.unique(selected_bins)) != len(selected_bins):
        raise ValueError('Distinct genomic bins required')
    selected_genes = np.flatnonzero(np.isin(gene_bins, selected_bins) & (fractions > 0))
    lookup = {int(b): j for j, b in enumerate(selected_bins)}
    local_bin = np.array([lookup[int(gene_bins[g])] for g in selected_genes], dtype=int)
    O = np.asarray(counts[:, selected_genes], float)
    E = np.asarray(expected_bins[:, local_bin], float)*fractions[selected_genes]
    r = np.broadcast_to(np.asarray(gene_sizes, float), (len(gene_bins),))[selected_genes]
    return O, E, r, selected_genes, local_bin


@contextmanager
def event_likelihood(registry, mode, trace, directory, *, region_provider=None):
    if mode not in ('control', 'moment_bin', 'gene'):
        raise ValueError(mode)
    original = clones._fit_event_mixture
    event_code = clones.SpecimenModel.event_mixtures.__code__
    cache = {}

    def fit(O, E, r, step, *, maxiter=1000):
        frame = sys._getframe(1)
        if frame.f_code is not event_code:
            return original(O, E, r, step, maxiter=maxiter)
        loc = frame.f_locals
        model, cells = loc['self'], np.asarray(loc['cells'], int)
        info = registry[id(model.rows)]
        if not model.thin or model.alleles is not None or model.graph is not None:
            raise ValueError('Gene proposals require gene-native thinning, expression-only and no spatial graph')
        if 'prop' not in info:
            raise RuntimeError('Proposal counts were not generated before event search')
        selected_bins = loc['cols'][loc['part']]
        source = f"event_chrom{loc['c']}_bins{selected_bins[0]}-{selected_bins[-1]}"
        key = (id(model.rows), cells.tobytes())
        if key not in cache:
            X = info['prop'][cells].toarray()
            sizes = info['resolved_gene_size']*info['share']
            fractions, bin_sizes, support = composition(X, info['gbin'], sizes, model.J)
            pool = directory/f'parent_{len(cache):04d}.npz'
            np.savez_compressed(pool, fractions=fractions, gene_sizes=sizes, bin_sizes=bin_sizes,
                bin_support=support, cells=cells, global_rows=model.rows.rows[cells], gbin=info['gbin'])
            cache[key] = dict(X=X, fractions=fractions, bin_sizes=bin_sizes, support=support,
                              sizes=sizes, path=pool)
        parent = cache[key]
        if region_provider is None:
            gene_O, gene_E, gene_r, genes, local_bin = region_arrays(parent['X'], info['gbin'],
                parent['fractions'], parent['sizes'], selected_bins, E)
        else:
            if mode != 'gene':
                raise ValueError('Reference gene composition requires the gene likelihood arm')
            gene_O, gene_E, gene_r, genes, local_bin = region_provider(model, info, cells, selected_bins, E)
        count_check, mean_check = np.zeros_like(O, dtype=float), np.zeros_like(E, dtype=float)
        for j in range(len(selected_bins)):
            count_check[:, j] = gene_O[:, local_bin == j].sum(axis=1)
            mean_check[:, j] = gene_E[:, local_bin == j].sum(axis=1)
        np.testing.assert_array_equal(count_check, O)
        supported = bool(parent['support'][selected_bins].all()) if region_provider is None else bool(
            np.isfinite(gene_E).all() and np.all(gene_E > 0))
        operations = 4*len(genes)+32
        unit = np.finfo(float).eps
        allowance = operations*unit/(1-operations*unit)*max(1.,float(np.max(np.abs(E))))
        error = float(np.max(np.abs(mean_check-E)))
        if supported and error > allowance:
            raise ArithmeticError('Gene expectations do not conserve native bin expectations')
        index = len(trace)
        payload = directory/f'event_{index:05d}.npz'
        actual = {} if region_provider is None else dict(gene_O=gene_O, gene_E=gene_E, gene_r=gene_r)
        np.savez_compressed(payload, native_O=O, native_E=E, native_r=r, selected_bins=selected_bins,
                            selected_genes=genes, local_bin=local_bin, **actual)
        if mode == 'control':
            result = original(O,E,r,step,maxiter=maxiter)
        elif not supported:
            model.log.append(dict(action='gene_event_unresolved', source=source,
                                  reason='no_discovery_counts_in_a_region_bin'))
            result = dict(gain=0., responsibility=np.full(len(O), .5),
                numerical_gain_budget=0., parameters=None, converged=False, iterations=0)
        elif mode == 'moment_bin':
            result = original(O,E,parent['bin_sizes'][selected_bins],step,maxiter=maxiter)
        else:
            result = original(gene_O,gene_E,gene_r,step,maxiter=maxiter)
        trace.append(dict(index=index, mode=mode, source=source, source_matrix=info['source_matrix'],
            case_seed=info['case_seed'], parent=str(parent['path']), arrays=str(payload),
            proposal_cells=len(cells), region_genes=len(genes), gene_composition_supported=supported,
            maximum_expectation_discrepancy=error, expectation_rounding_allowance=float(allowance),
            gain=float(result['gain']), numerical_gain_budget=float(result['numerical_gain_budget']),
            responsibility=result['responsibility'].tolist(),
            parameters=None if result['parameters'] is None else result['parameters'].tolist(),
            nuisance_source='proposal_counts_only' if region_provider is None else 'fit_reference_gene_composition',
            confirmation_counts_used=0))
        return result

    clones._fit_event_mixture = fit
    try:
        yield
    finally:
        clones._fit_event_mixture = original
