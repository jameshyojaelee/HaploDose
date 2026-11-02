"""End-to-end run: two streaming passes over the counts, then one fit per specimen.

Pass 1: per-gene totals (gene filter). Pass 2: per-cell bin sums written to a
disk-backed cells x bins matrix. Everything after reads that matrix in chunks.
Optional phased SNP counts (`alleles`) are read once into sparse cells x SNPs
matrices; each specimen uses the SNPs that pass its patient's reference cells.
"""
from __future__ import annotations

import platform
import resource
import sys
import time
from dataclasses import asdict
from importlib import metadata

import jax
import numpy as np
import pandas as pd
from scipy import sparse

from . import __version__
from .alleles import (AlleleParams, SpecimenAlleles, balanced_cells, joint_state, normal_snp_filter,
                      read_alleles, segment_alleles)
from .alleles import cell_loglik as allele_loglik
from ._native import bin_sums, grouped_column_sums, native_available
from .bins import RIBOSOMAL, arm_labels, arm_starts, chromosome_starts, make_bins
from .clones import CloneParams, Rows, clone_order, fit_specimen, neighbour_graph, reference_depth_slopes, tree_table
from .io import (CountSource, new_output_dir, open_memmap, read_cells, read_gene_list,
                 read_positions, read_site_list, write_json, write_table)
from .kernels.estep import cell_loglik, posteriors, programme_scores
from .profile import log_ratio, segments_table
from .programmes import learn_programmes

QUERY_ROLES = ('tumour_epithelium', 'query')


def run(counts, cells, positions, out, exclude_genes=None, bin_genes=10, min_gene_mean=0.02,
        chunk=4096, query_roles=QUERY_ROLES, specimens=None, params: CloneParams | None = None,
        exclude_ribosomal=True, alleles=None, allele_exclude_genes=None, allele_exclude_sites=None,
        allele_params: AlleleParams | None = None):
    params = params or CloneParams()
    allele_params = allele_params or AlleleParams()
    t0 = time.monotonic()
    started = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    out = new_output_dir(out)
    with CountSource(counts) as src:
        meta = read_cells(cells, src.cell_ids)
        pos = read_positions(positions)
        exclude = read_gene_list(exclude_genes)
        query = meta.role.isin(query_roles).to_numpy()
        ref = meta.role.eq('reference').to_numpy()
        fit = ref & meta.ref_split.eq('fit').to_numpy()
        hold = ref & meta.ref_split.eq('holdout').to_numpy()
        if fit.sum() < 20:
            raise ValueError('Need at least 20 fit reference cells')
        # Pass 1: gene sums per fit reference type, per dataset's query cells, and held-out references.
        fit_types = sorted(set((meta.dataset + '|' + meta.ref_group)[fit]))
        datasets = sorted(set(meta.dataset[query | ref]))
        tkey = (meta.dataset + '|' + meta.ref_group).to_numpy()
        group = np.full(len(meta), -1, np.int32)
        group[fit] = pd.Index(fit_types).get_indexer(tkey[fit])
        group[query] = len(fit_types) + pd.Index(datasets).get_indexer(meta.dataset[query])
        group[hold] = len(fit_types) + len(datasets)
        n_groups = len(fit_types) + len(datasets) + 1
        sums = np.zeros((n_groups, src.n_genes))
        for a, b, X in src.chunks(chunk):
            sums += grouped_column_sums(X, group[a:b], n_groups)
        gene_mean = sums.sum(axis=0) / max(1, int((query | ref).sum()))
        mismatch = type_mismatch_genes(sums, fit_types, datasets, src.genes)
        bins, gene_bin, used = make_bins(src.genes, pos, gene_mean, exclude | mismatch, bin_genes, min_gene_mean,
                                         exclude_pattern=RIBOSOMAL if exclude_ribosomal else None)
        J = len(bins)
        Y = open_memmap(out / 'bin_counts.npy', np.float32, (src.n_cells, J))
        library = np.empty(src.n_cells)
        for a, b, X in src.chunks(chunk):
            Y[a:b], library[a:b] = bin_sums(X, gene_bin, J)
        Y.flush()
        n_genes_input = src.n_genes
        # Per dataset: fit-reference counts of each used gene, to split a bin's expectation across its genes.
        ref_gene = {ds: sums[[t for t, key in enumerate(fit_types) if key.startswith(ds + '|')]].sum(axis=0)
                    [used.column.to_numpy()] for ds in datasets}
    Lbin = np.concatenate([Y[a:a + chunk].sum(axis=1, dtype=np.float64) for a in range(0, len(meta), chunk)])
    allele_data = None if alleles is None else read_alleles(
        alleles, meta.cell_id.to_numpy(), read_gene_list(allele_exclude_genes), read_site_list(allele_exclude_sites))
    reset = chromosome_starts(bins)
    write_table(bins, out / 'bins.tsv')
    write_table(used, out / 'genes.tsv.gz')
    write_table(pd.DataFrame({'gene': sorted(mismatch)}), out / 'type_mismatch_genes.tsv')

    ref_key = (meta.dataset + '|' + meta.ref_group).to_numpy()
    types = sorted(set(ref_key[fit & (Lbin > 0)]))
    mu, phi_parts = reference_rates(Y, Lbin, ref_key, fit, types, chunk, params.seed)

    todo = sorted(meta.specimen[query].unique()) if specimens is None else list(specimens)
    status = np.where(query, 'not_fitted', 'not_query').astype(object)
    status[query & (Lbin < params.min_umi)] = 'insufficient_depth'
    clone = np.full(len(meta), '', dtype=object)
    post_max = np.full(len(meta), np.nan)
    p_normal = np.full(len(meta), np.nan)
    margin = np.full(len(meta), np.nan)
    programmes, depth_gamma = {}, {}
    post_store, tables = {}, {k: [] for k in ['clones', 'segments', 'profiles', 'differences', 'log', 'tree',
                                             'reference_check', 'specimens', 'allele_segments', 'arms']}
    allele_depth = np.zeros(len(meta))
    in_sample_ref = np.zeros(len(meta), bool)
    allele_snps = np.zeros(len(meta), dtype=int)
    for spec in todo:
        rows = np.flatnonzero(query & meta.specimen.eq(spec).to_numpy() & (Lbin >= params.min_umi))
        if len(rows) < params.min_cells:
            status[query & meta.specimen.eq(spec).to_numpy() & (Lbin >= params.min_umi)] = 'too_few_cells'
            continue
        dataset = meta.dataset.iloc[rows[0]]
        cols = [t for t, key in enumerate(types) if key.startswith(dataset + '|')]
        if not cols:
            raise ValueError(f'No fit reference cells for dataset {dataset}')
        phi = phi_parts[0][:, cols].sum(axis=1) / np.maximum(phi_parts[1][:, cols].sum(axis=1), 1e-12)
        phi = np.clip(phi, 1e-3, 10.0)
        if dataset not in programmes:
            programmes[dataset] = None if params.n_programmes <= 0 else learn_programmes(
                Y, Lbin, mu[:, cols], phi,
                [np.flatnonzero(fit & (ref_key == types[c]) & (Lbin >= params.min_umi)) for c in cols],
                params.n_programmes, arm_starts(bins) if params.programme_window <= 0 else reset,
                params.seed, params.programme_window)
        prog = programmes[dataset]
        if params.depth_source == 'reference' and dataset not in depth_gamma:
            depth_gamma[dataset] = reference_depth_slopes(
                Y, Lbin, [np.flatnonzero(fit & (ref_key == types[c]) & (Lbin >= params.min_umi)) for c in cols],
                params.min_cells, chunk)
        spec_alleles, snp_cols = None, None
        same_patient = meta.patient.eq(meta.patient.iloc[rows[0]]).to_numpy()
        if allele_data is not None:
            keep, kappa, n_biased = normal_snp_filter(allele_data.A, allele_data.D, np.flatnonzero(fit & same_patient),
                                                      allele_params)
            covered = np.asarray(allele_data.D[rows].sum(axis=0)).ravel() > 0
            snp_cols = np.flatnonzero(keep & covered)
            spec_alleles = SpecimenAlleles(allele_data.A[rows][:, snp_cols], allele_data.D[rows][:, snp_cols],
                                           allele_data.snps.iloc[snp_cols], kappa, allele_params, n_biased)
            spec_alleles.bins = bins
            allele_depth[rows] = np.asarray(spec_alleles.D.sum(axis=1)).ravel()
            allele_snps[rows] = np.diff(spec_alleles.D.indptr)
        # In-sample normal reference: cells balanced at the specimen's strong allelic imbalance.
        mu_spec, type_keys, n_in_sample = mu[:, cols], [types[c] for c in cols], 0
        if spec_alleles is not None and params.in_sample_reference:
            cand, _ = balanced_cells(spec_alleles)
            if len(cand) >= params.min_cells:
                crow = rows[cand]
                total = np.zeros(J)
                for a in range(0, len(crow), chunk):
                    total += np.asarray(Y[crow[a:a + chunk]], np.float64).sum(axis=0)
                rate = np.maximum(total / Lbin[crow].sum(), 0.1 / Lbin[crow].sum())
                mu_spec = np.column_stack([mu_spec, rate / rate.sum()])
                type_keys = type_keys + [f'{dataset}|in_sample_normal']
                n_in_sample = len(cand)
                in_sample_ref[crow] = True
        graph = None
        if {'x', 'y'} <= set(meta.columns) and params.spatial_neighbors > 0:
            xy = meta[['x', 'y']].apply(pd.to_numeric, errors='coerce').to_numpy(np.float64)[rows]
            if np.isfinite(xy).all():
                graph = neighbour_graph(xy, params.spatial_neighbors)
        res = fit_specimen(Rows(Y, rows, chunk), Lbin[rows], mu_spec, reset, phi, params, prog, spec_alleles,
                           graph, gene_loader(counts, used, rows), depth_gamma.get(dataset))
        summarise(spec, rows, res, bins, reset, params, type_keys, tables, status, clone,
                  post_max, p_normal, margin, post_store,
                  genes=(gene_loader(counts, used, rows), gene_shares(ref_gene[dataset], used.bin.to_numpy(), J)))
        tables['specimens'][-1]['in_sample_reference_cells'] = n_in_sample
        hold_rows = np.flatnonzero(hold & meta.dataset.eq(dataset).to_numpy() & (Lbin >= params.min_umi))
        hold_alleles = None
        if spec_alleles is not None:
            # Only the patient's own held-out cells carry its heterozygous SNPs.
            A = allele_data.A[hold_rows][:, snp_cols].multiply(same_patient[hold_rows][:, None]).tocsr()
            D = allele_data.D[hold_rows][:, snp_cols].multiply(same_patient[hold_rows][:, None]).tocsr()
            hold_alleles = (pd.Index(hold_rows), A, D)
        tables['reference_check'].append(reference_check(spec, hold_rows, ref_key, types, cols, Y, Lbin, mu,
                                                         phi, res, params, prog, hold_alleles=hold_alleles))

    k_max = max([p.shape[1] for _, p in post_store.values()], default=1)
    P = open_memmap(out / 'posteriors.npy', np.float16, (len(meta), k_max), fill=np.nan)
    for rows, p in post_store.values():
        P[rows, :p.shape[1]] = p.astype(np.float16)
    P.flush()
    cells_out = meta.assign(library_UMI=library, binned_UMI=Lbin, allele_depth=allele_depth, allele_snps=allele_snps,
                            in_sample_reference=in_sample_ref,
                            status=status, clone=clone,
                            clone_posterior=post_max, p_normal=p_normal, assignment_margin=margin,
                            malignant_call=(query & (status == 'assigned')
                                            & ~pd.Series(clone).str.endswith('_c0').to_numpy()
                                            & (p_normal < 0.5)))
    write_table(cells_out, out / 'cells.tsv.gz')
    for name, rows in tables.items():
        if name == 'allele_segments' and allele_data is None:
            continue
        write_table(pd.DataFrame(rows), out / f'{name}.tsv.gz' if name == 'profiles' else out / f'{name}.tsv')
    qc = dict(haplodose=__version__, started=started, wall_seconds=time.monotonic() - t0,
              max_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2,
              n_cells=len(meta), n_genes_input=n_genes_input, n_genes_used=len(used), n_bins=J,
              n_query=int(query.sum()), n_fit_reference=int(fit.sum()), n_holdout_reference=int(hold.sum()),
              reference_groups=types, specimens_fitted=len(post_store),
              programme_variance_share={d: None if p is None else p.variance_share for d, p in programmes.items()},
              alleles=None if allele_data is None else dict(
                  path=str(alleles), n_snps=int(allele_data.A.shape[1]), excluded_rows=allele_data.excluded,
                  exclude_genes=str(allele_exclude_genes), exclude_sites=str(allele_exclude_sites),
                  params=asdict(allele_params)),
              native_kernels=native_available(),
              backend=jax.default_backend(), devices=[str(d) for d in jax.devices()],
              host=platform.node(), python=sys.version.split()[0], versions=versions(),
              params=dict(asdict(params), bin_genes=bin_genes, min_gene_mean=min_gene_mean, chunk=chunk,
                          query_roles=list(query_roles), exclude_genes=str(exclude_genes),
                          n_excluded_names=len(exclude), exclude_ribosomal=exclude_ribosomal),
              n_type_mismatch_genes=len(mismatch),
              note='Relative RNA dosage states and RNA-inferred clones; not absolute DNA copy number '
                   'or validated DNA clones.')
    write_json(qc, out / 'qc.json')
    return qc


def type_mismatch_genes(sums, fit_types, datasets, genes, max_log_ratio=np.log(5.0), min_counts=20):
    """Genes whose query expression no reference type matches within a dosage-plausible ratio.

    Copy number moves a gene by at most beta * log(cn / 2) (~0.8 at cn 8, beta 0.6). A gene
    whose median-centred query/reference log ratio exceeds log 5 against every reference
    type of its dataset reflects cell type, not dosage (e.g. epithelial genes when only
    immune and stromal references exist); such genes would dominate their bins."""
    out = set()
    names = pd.Index(np.asarray(genes, dtype=str)).str.upper()
    for d, ds in enumerate(datasets):
        q = sums[len(fit_types) + d]
        cols = [t for t, key in enumerate(fit_types) if key.startswith(ds + '|')]
        if not cols or q.sum() == 0:
            continue
        enough = q >= min_counts
        best = np.full(len(q), np.inf)
        for t in cols:
            rf = sums[t]
            r = np.log((q + 0.5) / q.sum()) - np.log((rf + 0.5) / rf.sum())
            r -= np.median(r[enough])
            best = np.minimum(best, np.abs(r))
        out |= set(names[enough & (best > max_log_ratio)])
    return out


def gene_shares(ref_counts, gbin, n_bins, prior=0.5):
    """Each used gene's share of its bin in the reference (shares sum to 1 per bin), and the gene bins."""
    w = np.asarray(ref_counts, np.float64) + prior
    tot = np.bincount(gbin, weights=w, minlength=n_bins)
    return w / tot[gbin], gbin


def gene_level_log_ratio(y, E_bin, share, gbin, r_bin, min_expected=1.0):
    """Gene-level bin log ratio (log_ratio_gene): median over the bin's genes of log((y_g + 0.5) / (E_g + 0.5)).

    E_g splits the bin's model expectation by reference gene shares. Summed-count bins let a few highly
    expressed genes dominate (Kinker single-line rho vs DepMap 0.852 summed vs 0.862 per-gene median,
    results/singleline_decomp_20261101_2354). Genes with E_g < min_expected are left out (their log ratio
    shrinks to 0 in small clones); bins with none keep the summed-count value. The result is shifted so
    its median over bins equals that of r_bin (the summed-count, centred profile). Reported next to
    log_ratio, not instead of it: U-CAN spatial rho 7/11 patients higher/lower (p 0.52) and implanted-region
    amplitude in clones < 100 cells 0.419 -> 0.411 (results/genelr_20261102_0001)."""
    E = E_bin[gbin] * share
    ok = E >= min_expected
    lr = pd.Series(np.log((y + 0.5) / (E + 0.5))[ok]).groupby(gbin[ok]).median()
    out = r_bin.copy()
    if len(lr):
        v = lr.to_numpy()
        out[lr.index.to_numpy()] = v - np.median(v) + np.median(r_bin[lr.index.to_numpy()])
    return out


def gene_loader(counts, used, rows, block=4096):
    """Loader of gene-level counts (used genes only) for specimen-local row positions.

    Rows are read in storage blocks of `block` rows and selected in memory: a random subsample has few
    consecutive rows, so reading only contiguous runs meant one small read per cell."""
    cols = used.column.to_numpy()
    gbin = used.bin.to_numpy()

    def load(local):
        glob = np.asarray(rows)[np.asarray(local)]
        order = np.argsort(glob, kind='stable')
        g = glob[order]
        parts = []
        with CountSource(counts) as src:
            for blk in np.unique(g // block):
                sel = g[g // block == blk]
                parts.append(src.read(int(sel[0]), int(sel[-1]) + 1)[sel - sel[0]][:, cols])
        out = sparse.vstack(parts).tocsr()
        return out[np.argsort(order, kind='stable')], gbin
    return load


def reference_rates(Y, Lbin, ref_key, fit, types, chunk, seed, max_cells=2000):
    """Bin rates per reference group (columns sum to 1) and NB moment parts for dispersion."""
    rng = np.random.default_rng(seed)
    J = Y.shape[1]
    mu = np.zeros((J, len(types)))
    num = np.zeros((J, len(types)))
    den = np.zeros((J, len(types)))
    for t, key in enumerate(types):
        rows = np.flatnonzero(fit & (ref_key == key) & (Lbin > 0))
        total = np.zeros(J)
        for a in range(0, len(rows), chunk):
            total += np.asarray(Y[rows[a:a + chunk]], np.float64).sum(axis=0)
        floor = 0.1 / max(Lbin[rows].sum(), 1.0)
        rate = np.maximum(total / max(Lbin[rows].sum(), 1.0), floor)
        mu[:, t] = rate / rate.sum()
        sub = np.sort(rng.choice(rows, min(len(rows), max_cells), replace=False))
        Ys = np.asarray(Y[sub], np.float64)
        m = Lbin[sub, None] * mu[None, :, t]
        num[:, t] = ((Ys - m) ** 2 - m).sum(axis=0)
        den[:, t] = (m ** 2).sum(axis=0)
    return mu, (num, den)


def summarise(spec, rows, res, bins, reset, params, type_keys, tables, status, clone, post_max, p_normal,
              margin, post_store, genes=None, chunk=4096):
    model = res['model']
    order = clone_order(res['R'])
    R, LL = res['R'][:, order], res['LL'][:, order]
    Yk, Lk, XLk = res['Yk'][order], res['Lk'][order], res['XLk'][order]
    prof = {k: v[order] for k, v in res['prof'].items()}
    names = [f'{spec}_c{i}' for i in range(len(order))]
    mapc = R.argmax(axis=1)
    status[rows] = 'assigned'
    clone[rows] = np.asarray(names, dtype=object)[mapc]
    post_max[rows] = R.max(axis=1)
    p_normal[rows] = R[:, 0]
    if LL.shape[1] > 1:
        top = np.sort(LL, axis=1)
        margin[rows] = top[:, -1] - top[:, -2]
    post_store[spec] = (rows, R)
    Yg = None
    if genes is not None:
        # Gene counts pooled per clone with the same responsibilities as the bin pools Yk.
        load, (share, gbin) = genes
        Yg = 0.0
        for a in range(0, len(rows), chunk):
            G, _ = load(np.arange(a, min(a + chunk, len(rows))))
            Yg = Yg + np.asarray(G.T @ R[a:a + chunk].astype(np.float64)).T
    chrom = bins.chromosome.to_numpy()
    nb = len(bins)
    al = model.alleles
    arm = arm_labels(bins)
    for k, name in enumerate(names):
        segs = segments_table(prof['post'][k], prof['path'][k], prof['dose'][k], Yk[k], prof['E'][k], reset,
                              model.beta, prof['tau'][k], model.pp)
        resolved = 0
        for s in segs:
            a, b = s['first_bin'], s['last_bin']
            resolved += s['resolved']
            row = dict(specimen=spec, clone=name, chromosome=chrom[a], start=int(bins.start.iloc[a]),
                       end=int(bins.end.iloc[b]), n_genes=int(bins.n_genes.iloc[a:b + 1].sum()), **s)
            if al is not None:
                summary = segment_alleles(model.allele_post[order[k]], al.snps, model.Dk[order[k]], bins, a, b)
                row.update(summary, joint_state=joint_state(s['cn_state'], summary['p_imbalanced'], summary['n_snps'],
                                                            summary['allele_depth'], al.params))
            tables['segments'].append(row)
        allele_resolved = 0
        if al is not None:
            for r in al.runs(model.allele_post[order[k]], model.Dk[order[k]]):
                allele_resolved += r['resolved']
                tables['allele_segments'].append(dict(specimen=spec, clone=name, **r))
        r = log_ratio(Yk[k], prof['E'][k]) - prof['centre'][k]
        # Extra column, existing log_ratio unchanged: per bin, median over genes of the gene log ratio.
        rg = r if Yg is None else gene_level_log_ratio(Yg[k], prof['E'][k], share, gbin, r)
        frame = pd.DataFrame({'arm': arm, 'r': r, 'rg': rg, 'dose': prof['dose'][k]})
        for a_name, g in frame.groupby('arm', sort=False):
            tables['arms'].append(dict(specimen=spec, clone=name, arm=a_name, n_bins=len(g),
                                       median_log_ratio=float(g.r.median()), median_log_ratio_gene=float(g.rg.median()),
                                       mean_log_dose=float(g.dose.mean())))
        tables['profiles'].extend(dict(specimen=spec, clone=name, bin=j, log_ratio=float(r[j]),
                                       log_ratio_gene=float(rg[j]),
                                       log_dose=float(prof['dose'][k][j]), cn_state=int(prof['path'][k][j] + 1))
                                  for j in range(nb))
        overlap = np.nan if k == 0 or model.programmes is None else model.programmes.overlap(prof['dose'][k],
                                                                                             model.base)
        tables['clones'].append(dict(specimen=spec, clone=name, column=k, normal_like=k == 0,
                                     n_cells_map=int((mapc == k).sum()), n_eff=float(R[:, k].sum()),
                                     pooled_UMI=float(Lk[k]), resolved_segments=resolved,
                                     tau=float(prof['tau'][k]), centre=float(prof['centre'][k]),
                                     programme_overlap=overlap,
                                     programme_shift=float(prof['programme_shift'][k]),
                                     allele_segments_resolved=allele_resolved if al is not None else np.nan))
    # Pairwise differences between clones with enough cells, on the held-out (test) bins, depth-matched.
    # The pools are MAP cells, so a clone must also have min_cells at MAP (a clone sharing most of its dose
    # with another can keep soft weight >= min_cells yet win no cell, which left an empty pool).
    ok = [k for k in range(len(names))
          if R[:, k].sum() >= params.min_cells and (mapc == k).sum() >= params.min_cells]
    dist = np.zeros((len(ok), len(ok)))
    pairs = [(i, j) for i in range(len(ok)) for j in range(i + 1, len(ok))]
    if pairs:
        A, B = [ok[i] for i, _ in pairs], [ok[j] for _, j in pairs]
        diff = model.matched_differences(mapc, list(zip(A, B)))   # depth-matched pooling of MAP cells
        for (i, j), d in zip(pairs, diff):
            genes = 0
            # Report only segments that pass the same gene and allele vetoes as merging and splitting
            # (vetoed segments are written to the log table by model.supported).
            a_, b_ = order[ok[i]], order[ok[j]]
            ad = (None,) * 4 if model.Ak is None else (model.Ak[a_], model.Dk[a_], model.Ak[b_], model.Dk[b_])
            for s in model.supported(d['segments'], *ad, labels=mapc, a=ok[i], b=ok[j]):
                a, b = s['first_bin'], s['last_bin']
                n = int(bins.n_genes.iloc[a:b + 1].sum())
                genes += n
                tables['differences'].append(dict(specimen=spec, clone_a=names[ok[i]], clone_b=names[ok[j]],
                                                  chromosome=chrom[a], start=int(bins.start.iloc[a]),
                                                  end=int(bins.end.iloc[b]), n_genes=n, **s))
            dist[i, j] = dist[j, i] = genes
    tables['tree'].extend(dict(specimen=spec, **r) for r in tree_table(dist, [names[k] for k in ok]))
    tables['log'].extend(dict(specimen=spec, **e) for e in model.log)
    tables['specimens'].append(dict(specimen=spec, n_cells=len(rows), n_clones=len(names) - 1,
                                    normal_like_cells=int((mapc == 0).sum()), beta=model.beta, tau=model.tau,
                                    centre=model.centre, median_phi=float(np.median(model.phi)),
                                    mixture=';'.join(f'{t}={w:.3f}' for t, w in zip(type_keys, model.weights)),
                                    n_snps=np.nan if al is None else len(al.snps),
                                    snps_biased_in_normals=np.nan if al is None else al.n_biased,
                                    allele_kappa=np.nan if al is None else al.kappa,
                                    spatial_graph=model.graph is not None, potts_lambda=model.lam,
                                    phase_switch_nu=np.nan if al is None else al.nu))


def reference_check(spec, hold_rows, ref_key, types, cols, Y, Lbin, mu, phi, res, params, prog=None,
                    per_type=300, hold_alleles=None):
    """Held-out normal cells scored on their own group baseline: clone 0 vs this specimen's clones.

    Tests the dosage evidence alone: the baseline is correct by construction and the
    prior is flat (0.5 copy-neutral). `expected` re-draws each tested cell from its own
    baseline under the fitted negative binomial and scores it the same way: the rate a
    correctly specified model would give. With alleles, the patient's own held-out
    cells add their SNP evidence (re-drawn as balanced for `expected`). Observed well above expected means the
    likelihood is overconfident; both near 0.5 means the profiles are too weak to
    separate single normal cells. With programmes, cells are scored as in the fit:
    scores around the first-pass posterior-weighted rate, then re-scored.
    """
    rng = np.random.default_rng(params.seed)
    doses = res['doses']
    # Flat prior: half on copy-neutral, half spread over this specimen's clones.
    pi = np.full(len(doses), 0.5 / max(1, len(doses) - 1))
    pi[0] = 0.5 if len(doses) > 1 else 1.0
    called = tested = expected = 0
    for c in cols:
        rows = hold_rows[ref_key[hold_rows] == types[c]]
        if not len(rows):
            continue
        rows = np.sort(rng.choice(rows, min(len(rows), per_type), replace=False))
        lr = np.log(mu[:, c])[None] + doses
        lr = lr - np.log(np.exp(lr).sum(axis=1, keepdims=True))
        allele = draw_allele = None
        model = res['model']
        if hold_alleles is not None and model.B is not None:
            index, A, D = hold_alleles
            pos = index.get_indexer(rows)
            Ar, Dr = A[pos], D[pos]
            allele = (Ar, Dr, model.B, model.alleles.params.cell_shrink)
            Ad = Dr.copy()
            Ad.data = rng.binomial(Dr.data.astype(np.int64), 0.5).astype(np.float32)
            draw_allele = (Ad, Dr, model.B, model.alleles.params.cell_shrink)
        P = _score(np.asarray(Y[rows], np.float32), Lbin[rows], lr, phi, pi, prog, allele)
        called += int((P[:, 0] < 0.5).sum())
        m = Lbin[rows, None] * mu[None, :, c]
        draw = rng.poisson(rng.gamma(1.0 / phi[None], m * phi[None])).astype(np.float32)
        P = _score(draw, draw.sum(axis=1), lr, phi, pi, prog, draw_allele)
        expected += int((P[:, 0] < 0.5).sum())
        tested += len(rows)
    return dict(specimen=spec, n_tested=tested, n_called_aneuploid=called,
                fraction_called_aneuploid=called / tested if tested else np.nan,
                expected_fraction=expected / tested if tested else np.nan)


def _score(Yc, Lc, lr, phi, pi, prog, allele=None):
    extra = 0.0 if allele is None else allele_loglik(*allele)
    P = posteriors(cell_loglik(Yc, Lc, lr, phi) + extra, np.log(pi))
    if prog is None:
        return P
    U = programme_scores(Yc, Lc, P.astype(np.float32), lr, phi, prog.W, prog.prior_prec)
    return posteriors(cell_loglik(Yc, Lc, lr, phi, U=U, W=prog.W) + extra, np.log(pi))


def versions():
    out = {}
    for pkg in ['jax', 'jaxlib', 'numpy', 'scipy', 'pandas', 'h5py', 'scanpy', 'anndata', 'leidenalg']:
        try:
            out[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            out[pkg] = None
    return out
