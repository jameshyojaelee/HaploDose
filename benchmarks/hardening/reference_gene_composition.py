"""Reference-only positive within-bin gene means, with declared depth transport."""
import sys

import numpy as np
from haplodose import pipeline


def table(counts, depth, groups, gene_bins, *, strata=5):
    counts, depth, groups, gene_bins = map(np.asarray,(counts,depth,groups,gene_bins))
    if (counts.ndim != 2 or counts.shape != (len(depth),len(gene_bins)) or groups.shape != depth.shape
            or not np.isfinite(counts).all() or np.any(counts < 0)
            or not np.isfinite(depth).all() or np.any(depth <= 0) or strata not in (1,5)):
        raise ValueError('Aligned fit-reference counts, positive observed depths and strata1/5 required')
    types = sorted(set(groups))
    nbin = int(gene_bins.max())+1
    if not np.array_equal(np.unique(gene_bins),np.arange(nbin)):
        raise ValueError('Every bin needs reference genes')
    fractions, knots, assignments = [], [], []
    for label in types:
        selected = np.flatnonzero(groups == label)
        x, exposure = np.log(depth[selected]), depth[selected]
        edges = np.quantile(x,np.linspace(0,1,strata+1)[1:-1])
        split = np.searchsorted(edges,x)
        per_type, centres = [], []
        for k in range(strata):
            subset = selected[split == k]
            if len(subset) < 2:
                raise ValueError('Reference composition stratum has fewer than two cells')
            weight = float(depth[subset].sum()/exposure.sum())
            gene = counts[subset].sum(axis=0)+.5*weight
            totals = np.bincount(gene_bins,weights=gene,minlength=nbin)
            per_type.append(gene/totals[gene_bins])
            centres.append(float(np.sum(np.log(depth[subset])*depth[subset])/depth[subset].sum()))
            assignments.append(dict(type=str(label),stratum=k,rows=subset.tolist(),regularizer_fraction=weight))
        if len(centres)>1 and not np.all(np.diff(centres)>0):
            raise ValueError('Distinct reference depth knots required')
        fractions.append(per_type)
        knots.append(centres)
    return dict(types=types,fraction=np.asarray(fractions),knots=np.asarray(knots),
        strata=strata,assignments=assignments)


def predict(fitted, depth, gene_bins, mixture_contributions):
    depth, gene_bins, mixture_contributions = np.asarray(depth,float),np.asarray(gene_bins,int),np.asarray(mixture_contributions,float)
    nbin = int(gene_bins.max())+1
    if mixture_contributions.shape != (len(fitted['types']),nbin):
        raise ValueError('Aligned native type contributions required')
    if (np.any(mixture_contributions < 0) or not np.isfinite(mixture_contributions).all()
            or np.any(mixture_contributions.sum(axis=0) <= 0) or np.any(depth <= 0)):
        raise ValueError('Supported type mixture and positive depths required')
    contribution = mixture_contributions/mixture_contributions.sum(axis=0)
    result = np.zeros((len(depth),len(gene_bins)))
    for t in range(len(fitted['types'])):
        knots, values = fitted['knots'][t],fitted['fraction'][t]
        if len(knots) == 1:
            mean = np.broadcast_to(values[0],result.shape)
        else:
            x = np.log(depth)
            lo = np.clip(np.searchsorted(knots,x)-1,0,len(knots)-2)
            weight = np.clip((x-knots[lo])/(knots[lo+1]-knots[lo]),0.,1.)
            mean = values[lo]*(1-weight[:,None])+values[lo+1]*weight[:,None]
        result += mean*contribution[t,gene_bins]
    sums = np.stack([result[:,gene_bins == b].sum(axis=1) for b in range(nbin)],axis=1)
    return result/sums[:,gene_bins]


def estimator(noise, mode, records, directory):
    original_factory = noise.estimator
    cache = {}
    def factory(noise_records, depth_strata=1):
        original = original_factory(noise_records,depth_strata)
        def size(info):
            value = original(info)
            if mode == 'control':
                return value
            frame = sys._getframe(1)
            while frame is not None and frame.f_code is not pipeline._run_single.__code__:
                frame = frame.f_back
            if frame is None:
                raise RuntimeError('Missing exact pipeline reference context')
            loc = frame.f_locals
            key = (str(loc['counts']),tuple(loc['used'].column),str(loc['dataset']))
            if key not in cache:
                meta = loc['meta']
                rows = np.flatnonzero(loc['fit'] & meta.dataset.eq(loc['dataset']).to_numpy())
                assert meta.iloc[rows].role.eq('reference').all()
                assert meta.iloc[rows].ref_split.eq('fit').all()
                X, gbin = pipeline.gene_loader(loc['counts'],loc['used'],rows)(np.arange(len(rows)))
                X = X.toarray().astype(float)
                assert np.array_equal(gbin,info['gbin'])
                depth = np.asarray(loc['Lbin'][rows],float)
                np.testing.assert_array_equal(X.sum(axis=1),depth)
                groups = loc['ref_key'][rows]
                fitted = table(X,depth,groups,gbin,strata=1 if mode == 'pooled' else 5)
                types = [loc['types'][c] for c in loc['cols']]
                assert fitted['types'] == types
                path = directory/f'reference_{len(cache):03d}.npz'
                mu = loc['mu'][:,loc['cols']].copy()
                np.savez_compressed(path,counts=X,depth=depth,groups=groups.astype(str),
                    global_rows=rows,cell_ids=meta.cell_id.to_numpy()[rows].astype(str),
                    gene_bins=gbin,mu=mu,fraction=fitted['fraction'],knots=fitted['knots'])
                cache[key] = dict(fitted=fitted,mu=mu,path=path)
                records.append(dict(source_counts=str(loc['counts']),reference=str(path),
                    reference_cells=len(rows),genes=len(gbin),mode=mode,
                    training_cell_ids=meta.cell_id.to_numpy()[rows].tolist(),
                    types=types,strata=fitted['assignments'],query_rows_used=0,holdout_rows_used=0))
            info['research_reference_composition'] = cache[key]
            return value
        return size
    noise.estimator = factory


def region_provider(records, directory):
    cache = {}
    def arrays(model, info, cells, selected_bins, expected_bins):
        reference = info['research_reference_composition']
        key = id(model.rows)
        if key not in cache:
            contribution = reference['mu'].T*model.weights[:,None]
            fractions = predict(reference['fitted'],model.L,info['gbin'],contribution)
            cache[key] = fractions
            path = directory/f'cell_composition_{len(cache)-1:03d}.npz'
            np.savez_compressed(path,fractions=fractions,depth=model.L,global_rows=model.rows.rows,
                contribution=contribution,gene_bins=info['gbin'])
            x, knots = np.log(model.L),reference['fitted']['knots']
            records.append(dict(cell_composition=str(path),reference=str(reference['path']),
                case_seed=info['case_seed'],cells=len(model.L),
                below_reference_knots=[int(np.sum(x < k[0])) for k in knots],
                above_reference_knots=[int(np.sum(x > k[-1])) for k in knots],
                tail_rule='hold_endpoint_composition; transport beyond support unvalidated'))
        fraction = cache[key]
        genes = np.flatnonzero(np.isin(info['gbin'],selected_bins))
        lookup = {int(b):j for j,b in enumerate(selected_bins)}
        local = np.array([lookup[int(info['gbin'][g])] for g in genes],int)
        O = info['prop'][cells][:,genes].toarray()
        E = expected_bins[:,local]*fraction[np.ix_(cells,genes)]
        r = info['resolved_gene_size'][genes]*info['share']
        return O,E,r,genes,local
    return arrays
