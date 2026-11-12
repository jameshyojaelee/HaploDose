"""Bounded known-noise gene thinning for the original panel; never a default."""
from contextlib import contextmanager
import inspect
from pathlib import Path

import numpy as np
from scipy import sparse
from haplodose import clones, pipeline


def csr_arrays(prefix, matrix):
    matrix = sparse.csr_matrix(matrix)
    return {prefix+'_data': matrix.data, prefix+'_indices': matrix.indices,
        prefix+'_indptr': matrix.indptr, prefix+'_shape': np.array(matrix.shape)}


@contextmanager
def gene_native(trace, directory, gene_size=10., *, case_seeds=(20261031, 20261107, 20261108, 20261109, 20261110)):
    original_fit, original_fold, original_gene = pipeline.fit_specimen, clones.Rows._fold, clones.thin_gene_counts
    signature = inspect.signature(original_fit)
    registry, current = {}, []

    def split(info):
        rows = info['rows']
        if rows.thin is None:
            raise RuntimeError('Gene fold requested before thinning setup')
        seed, _, share = rows.thin
        identity = (seed, share)
        if 'split_identity' in info:
            if identity != info['split_identity']:
                raise RuntimeError('Cannot change a cached research split')
            return
        gene_counts = info['X'].toarray()
        resolved_size = gene_size(info) if callable(gene_size) else gene_size
        sizes = np.broadcast_to(np.asarray(resolved_size, float), (gene_counts.shape[1],)).copy()
        if not np.isfinite(sizes).all() or np.any(sizes <= 0):
            raise ValueError('Gene thinning requires finite positive estimated sizes')
        info['resolved_gene_size'] = sizes
        prop = clones.thin_counts(gene_counts, rows.rows, sizes, seed, share)
        assert np.all(prop >= 0) and np.all(prop <= gene_counts)
        info.update(split_identity=identity, prop=sparse.csr_matrix(prop),
            prop_bins=np.asarray(prop@info['aggregation']), seed=seed, share=share)
        assert np.array_equal(np.asarray(gene_counts@info['aggregation']), info['full_bins'])

    def fold(rows, Yc, keys, name):
        if name is None:
            return original_fold(rows, Yc, keys, name)
        if id(rows) not in registry:
            raise RuntimeError('Unregistered gene-thinning rows; no bin fallback allowed')
        info = registry[id(rows)]
        split(info)
        local = np.array([info['lookup'][int(key)] for key in keys])
        assert np.array_equal(Yc, info['full_bins'][local])
        info['fold_calls'] += 1
        if name == 'prop':
            return info['prop_bins'][local].astype(np.float32)
        if name == 'test':
            return (Yc-info['prop_bins'][local]).astype(np.float32)
        raise ValueError('Unsupported research fold')

    def genes(Xg, gbin, Ybin, Ytest, seed):
        if not current:
            raise RuntimeError('Gene veto thinning outside active specimen')
        info = current[-1]
        split(info)
        local = info['loaded_local']
        assert seed == info['seed'] and np.array_equal(gbin, info['gbin'])
        assert (sparse.csr_matrix(Xg)-info['X'][local]).nnz == 0
        assert np.array_equal(Ybin, info['full_bins'][local])
        result = (info['X'][local]-info['prop'][local]).astype(np.float64).tocsr()
        assert np.array_equal((result@info['aggregation']).toarray(), Ytest)
        info['veto_local'], info['veto_genes'] = local.copy(), result
        info['veto_calls'] += 1
        return result

    def fit(*args, **kwargs):
        bound = signature.bind(*args, **kwargs)
        bound.apply_defaults()
        arguments = bound.arguments
        rows, params, loader = arguments['rows'], arguments['params'], arguments['gene_loader']
        if params.split_mode != 'thin' or params.split_null_draws or loader is None:
            raise ValueError('Known-gene prototype requires thin, gene evidence, and no unsupported split null')
        if arguments['alleles'] is not None or arguments['graph'] is not None:
            raise ValueError('Known-gene prototype is expression-only and nonspatial')
        X, gbin = loader(np.arange(rows.n))
        X, gbin = sparse.csr_matrix(X, dtype=np.float64), np.asarray(gbin, int)
        J = rows.Y.shape[1]
        assert np.all((gbin >= 0) & (gbin < J))
        aggregation = sparse.csr_matrix((np.ones(len(gbin)), (np.arange(len(gbin)), gbin)), shape=(len(gbin), J))
        full = (X@aggregation).toarray()
        assert np.array_equal(full, rows.take(np.arange(rows.n)))
        matrix_path = Path(rows.Y.filename)
        cases = [int(part.removeprefix('seed_')) for part in matrix_path.parts if part.startswith('seed_')]
        assert len(cases) == 1 and cases[0] in case_seeds
        info = dict(index=len(registry), rows=rows, X=X, gbin=gbin, aggregation=aggregation,
            case_seed=cases[0], source_matrix=str(matrix_path),
            full_bins=full, lookup={int(key): i for i, key in enumerate(rows.rows)},
            fold_calls=0, veto_calls=0)
        registry[id(rows)] = info

        def load(local):
            info['loaded_local'] = np.asarray(local, int).copy()
            return loader(local)

        arguments['gene_loader'] = load
        current.append(info)
        try:
            return original_fit(*bound.args, **bound.kwargs)
        finally:
            assert current.pop() is info

    pipeline.fit_specimen, clones.Rows._fold, clones.thin_gene_counts = fit, fold, genes
    try:
        yield registry
    finally:
        pipeline.fit_specimen, clones.Rows._fold, clones.thin_gene_counts = original_fit, original_fold, original_gene
        for info in registry.values():
            if 'prop' not in info:
                continue
            path = directory/f"model_{info['index']}.npz"
            arrays = dict(global_rows=info['rows'].rows, gbin=info['gbin'], full_bins=info['full_bins'],
                prop_bins=info['prop_bins'], seed=np.array(info['seed']), share=np.array(info['share']),
                gene_size=info['resolved_gene_size'] if callable(gene_size) else np.array(gene_size),
                **csr_arrays('source', info['X']), **csr_arrays('prop', info['prop']))
            if 'veto_genes' in info:
                arrays.update(veto_local=info['veto_local'], **csr_arrays('veto', info['veto_genes']))
            np.savez_compressed(path, **arrays)
            trace.append(dict(model=info['index'], cells=info['rows'].n, genes=info['X'].shape[1],
                case_seed=info['case_seed'], source_matrix=info['source_matrix'],
                bins=info['full_bins'].shape[1], seed=info['seed'], share=info['share'],
                oracle_gene_size=None if callable(gene_size) else gene_size,
                fold_calls=info['fold_calls'], veto_calls=info['veto_calls'], path=str(path)))
