"""End-to-end run: two streaming passes over the counts, then one fit per specimen.

Pass 1: per-gene totals (gene filter). Pass 2: per-cell bin sums written to a
disk-backed cells x bins matrix. Everything after reads that matrix in chunks.
Optional phased SNP counts (`alleles`) are read once into sparse cells x SNPs
matrices; each specimen uses the SNPs that pass its patient's reference cells.
"""
from __future__ import annotations

import platform
import json
import hashlib
import resource
import sys
import time
from dataclasses import asdict
from importlib import metadata
from pathlib import Path

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
from .device import select as select_device
from .clones import CloneParams, Rows, clone_order, fit_specimen, neighbour_graph, reference_depth_slopes, tree_table
from .io import (CountSource, new_output_dir, open_memmap, read_cells, read_gene_list,
                 read_positions, read_site_list, write_json, write_table, validate_input_contract)
from .provenance import start_run, finish_run, file_record
from .confirmation import Contrast, FrozenDiscovery, assign_roles, score_confirmation, contrast_diagnostic
from .datasets import materialize_datasets
from .kernels.estep import cell_loglik, posteriors, programme_scores
from .profile import log_ratio, segments_table
from .programmes import learn_programmes

QUERY_ROLES = ('tumour_epithelium', 'query')

# An empty result is still a typed table, rather than a headerless file that downstream readers
# cannot distinguish from a truncated output. Optional diagnostics may add columns to populated tables.
EMPTY_TABLE_COLUMNS = {
    'clones': ['specimen', 'clone', 'column', 'normal_like', 'n_cells_map', 'n_eff', 'pooled_UMI', 'pooled_count', 'count_units',
               'resolved_segments', 'tau', 'centre', 'programme_overlap', 'programme_shift', 'allele_segments_resolved'],
    'segments': ['specimen', 'clone', 'chromosome', 'start', 'end', 'n_genes', 'first_bin', 'last_bin',
                 'n_bins', 'cn_state', 'state', 'log_dose', 'rna_log_effect', 'relative_state_log_ratio',
                 'observed_mean_rna_log_ratio', 'rna_effect_estimator',
                 'posterior', 'pooled_UMI', 'expected_UMI', 'pooled_count', 'expected_count', 'count_units',
                 'se_log_ratio', 'se_rna_log_ratio_approx',
                 'three_sigma_rna_log_ratio_approx', 'uncertainty_method', 'uncertainty_scope',
                 'uncertainty_estimand', 'uncertainty_statistic', 'uncertainty_center_field',
                 'uncertainty_status', 'posterior_effect_interval_validated',
                 'uncertainty_calibrated', 'uncertainty_scale_z', 'scaled_se_rna_log_ratio_approx', 'resolved'],
    'profiles': ['specimen', 'clone', 'bin', 'log_ratio', 'log_ratio_gene', 'log_dose', 'assignment_dose', 'cn_state'],
    'differences': ['specimen', 'clone_a', 'clone_b', 'chromosome', 'start', 'end', 'n_genes', 'first_bin',
                    'last_bin', 'n_bins', 'direction', 'mean_diff', 'z', 'posterior',
                    'min_bin_posterior', 'posterior_rule'],
    'log': ['specimen', 'action'],
    'tree': ['specimen', 'node', 'left', 'right', 'genes_differing', 'n_clones'],
    'reference_check': ['specimen', 'n_tested', 'n_called_aneuploid', 'fraction_called_aneuploid', 'expected_fraction'],
    'specimens': ['specimen', 'n_cells', 'n_clones', 'normal_like_cells', 'beta', 'tau', 'centre', 'median_phi',
                  'mixture', 'n_snps', 'snps_biased_in_normals', 'allele_kappa', 'spatial_graph', 'potts_lambda',
                  'phase_switch_nu', 'in_sample_reference_cells'],
    'allele_segments': ['specimen', 'clone', 'chromosome', 'start', 'end', 'first_snp', 'last_snp', 'n_snps',
                        'depth', 'p_imbalanced', 'theta', 'resolved'],
    'arms': ['specimen', 'clone', 'arm', 'n_bins', 'median_log_ratio', 'median_log_ratio_gene', 'mean_log_dose'],
}


def check_identities(meta, query):
    """Each fitted specimen must belong to one patient and one dataset (patient picks the allele reference and
    dataset the expression reference); query cells need nonempty identifiers."""
    for col in ['specimen', 'patient', 'dataset']:
        empty = query & (meta[col].str.strip() == '').to_numpy()
        if empty.any():
            raise ValueError(f'{int(empty.sum())} query cells have an empty {col} '
                             f'(first: {meta.cell_id[empty].iloc[0]})')
    n = meta[query].groupby('specimen')[['patient', 'dataset']].nunique()
    mixed = n[(n > 1).any(axis=1)]
    if len(mixed):
        raise ValueError(f'specimens spanning more than one patient or dataset: {list(mixed.index[:5])}; give '
                         'each patient/dataset its own specimen name')


def run(counts, cells, positions, out, exclude_genes=None, bin_genes=10, min_gene_mean=0.02,
        chunk=4096, query_roles=QUERY_ROLES, specimens=None, params: CloneParams | None = None,
        exclude_ribosomal=True, alleles=None, allele_exclude_genes=None, allele_exclude_sites=None,
        allele_params: AlleleParams | None = None, device='cpu', *, counts_kind=None, genome_build=None,
        coordinate_system=None, gene_namespace=None, centromeres=None, max_gap_bp=None, discovery_fraction=None,
        allele_genome_build=None):
    """Schema-2 execution with declared measurement units and independent datasets.

    Each dataset has its own feature universe and fit. The parent manifest indexes
    subruns when more than one dataset is supplied. All groups remain exploratory.
    """
    contract = validate_input_contract(counts_kind, genome_build, coordinate_system, gene_namespace, centromeres)
    if alleles is not None:
        if allele_genome_build is None:
            raise ValueError('allele_genome_build must be explicitly declared for allele inputs')
        if allele_genome_build != genome_build:
            raise ValueError('Allele genome build must match the expression annotation genome build')
        contract['allele_genome_build'] = allele_genome_build
        contract['allele_coordinate_system'] = 'one-based-SNP-POS'
    if chunk <= 0 or bin_genes <= 0 or min_gene_mean < 0:
        raise ValueError('chunk and bin_genes must be positive; min_gene_mean must be nonnegative')
    with CountSource(counts, gene_namespace=gene_namespace) as source:
        meta = read_cells(cells, source.cell_ids)
    check_identities(meta, meta.role.isin(query_roles).to_numpy())
    params, allele_params = params or CloneParams(), allele_params or AlleleParams()
    if discovery_fraction is not None and (not np.isfinite(discovery_fraction) or not 0 < discovery_fraction < 1):
        raise ValueError('discovery_fraction must be strictly between 0 and 1')
    if discovery_fraction is not None and 'analysis_split' in meta:
        raise ValueError('Use explicit analysis_split metadata or discovery_fraction, not both')
    if alleles is not None and params.split_null_draws > 0:
        raise ValueError('allele-enabled split null is unsupported: observed and null statistics must match')
    if params.split_null_draws > 0:
        raise ValueError('Gene-veto-enabled split null is unsupported: observed and simulated searches must '
                         'apply the same gene evidence. Bin-only diagnostics remain available separately.')
    arguments = dict(exclude_genes=exclude_genes, bin_genes=bin_genes, min_gene_mean=min_gene_mean,
                     chunk=chunk, query_roles=query_roles, specimens=specimens, params=params,
                     exclude_ribosomal=exclude_ribosomal, alleles=alleles, allele_exclude_genes=allele_exclude_genes,
                     allele_exclude_sites=allele_exclude_sites, allele_params=allele_params, device=device,
                     counts_kind=counts_kind, genome_build=genome_build, coordinate_system=coordinate_system,
                     gene_namespace=gene_namespace, centromeres=centromeres, max_gap_bp=max_gap_bp,
                     discovery_fraction=discovery_fraction, allele_genome_build=allele_genome_build)
    out = new_output_dir(out)
    record = start_run(out, dict(counts=counts, cells=cells, positions=positions, alleles=alleles,
                                exclude_genes=exclude_genes, allele_exclude_genes=allele_exclude_genes,
                                allele_exclude_sites=allele_exclude_sites, centromeres=centromeres),
                       {**asdict(params), 'alleles': asdict(allele_params), 'bin_genes': bin_genes,
                        'min_gene_mean': min_gene_mean, 'chunk': chunk, 'device': device,
                        'query_roles': list(query_roles), 'specimens': specimens,
                        'exclude_ribosomal': exclude_ribosomal,
                        'max_gap_bp': max_gap_bp, 'discovery_fraction': discovery_fraction}, contract, mode='single_cell')
    try:
        if meta.dataset.nunique() == 1:
            qc = _run_single(counts, cells, positions, out, **arguments)
        else:
            children = materialize_datasets(counts, meta, out / 'dataset_inputs',
                                            gene_namespace=gene_namespace, chunk=chunk)
            summaries = []
            for child in children:
                child_out = out / child['key']
                summary = run(child['counts'], child['cells'], positions, child_out, **arguments)
                summaries.append(summary)
                child['output'] = str(child_out)
                child['manifest'] = str(child_out / 'run_manifest.json')
            write_json([{k: str(v) if isinstance(v, Path) else v for k, v in c.items()} for c in children],
                       out / 'datasets.json')
            record['children'] = [json.loads(Path(c['manifest']).read_text())['output']['id'] for c in children]
            qc = {k: sum(s[k] for s in summaries) for k in
                  ['n_cells', 'n_bins', 'n_query', 'n_fit_reference', 'n_holdout_reference', 'specimens_fitted',
                   'n_discovery_query', 'n_confirmation_query', 'n_confirmation_scored', 'n_confirmation_unresolved']}
            qc.update(schema_version=2, layout='independent_dataset_subruns', n_genes_used=None,
                      wall_seconds=sum(s['wall_seconds'] for s in summaries),
                      max_rss_gib=max(s['max_rss_gib'] for s in summaries), backend=summaries[0]['backend'])
            write_json(qc, out / 'qc.json')
        finish_run(out, record, backend=qc['backend'])
        return qc
    except Exception as error:
        finish_run(out, record, error=error)
        raise


def _run_single(counts, cells, positions, out, exclude_genes=None, bin_genes=10, min_gene_mean=0.02,
        chunk=4096, query_roles=QUERY_ROLES, specimens=None, params: CloneParams | None = None,
        exclude_ribosomal=True, alleles=None, allele_exclude_genes=None, allele_exclude_sites=None,
        allele_params: AlleleParams | None = None, device='cpu', *, counts_kind=None, genome_build=None,
        coordinate_system=None, gene_namespace=None, centromeres=None, max_gap_bp=None, discovery_fraction=None,
        allele_genome_build=None):
    """device: 'cpu' (default), 'gpu' or 'auto' for the JAX kernels (see haplodose.device)."""
    backend = select_device(device)
    params = params or CloneParams()
    allele_params = allele_params or AlleleParams()
    t0 = time.monotonic()
    started = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    out = Path(out)
    with CountSource(counts, gene_namespace=gene_namespace) as src:
        meta = read_cells(cells, src.cell_ids)
        pos = read_positions(positions, genome_build=genome_build, coordinate_system=coordinate_system,
                             gene_namespace=gene_namespace, centromeres=centromeres)
        exclude = read_gene_list(exclude_genes, gene_namespace)
        query = meta.role.isin(query_roles).to_numpy()
        meta = analysis_roles(meta, query, params, discovery_fraction)
        discovery = query & meta.analysis_split.eq('discovery').to_numpy()
        confirmation = query & meta.analysis_split.eq('confirmation').to_numpy()
        ref = meta.role.eq('reference').to_numpy()
        fit = ref & meta.ref_split.eq('fit').to_numpy()
        hold = ref & meta.ref_split.eq('holdout').to_numpy()
        check_identities(meta, query)
        if fit.sum() < 20:
            raise ValueError('Need at least 20 fit reference cells')
        # Pass 1: gene sums per fit reference type, per dataset's query cells, and held-out references.
        fit_types = sorted(set((meta.dataset + '|' + meta.ref_group)[fit]))
        datasets = sorted(set(meta.dataset[query | ref]))
        tkey = (meta.dataset + '|' + meta.ref_group).to_numpy()
        group = np.full(len(meta), -1, np.int32)
        group[fit] = pd.Index(fit_types).get_indexer(tkey[fit])
        group[discovery] = len(fit_types) + pd.Index(datasets).get_indexer(meta.dataset[discovery])
        n_groups = len(fit_types) + len(datasets) + 1
        sums = np.zeros((n_groups, src.n_genes))
        for a, b, X in src.chunks(chunk):
            sums += grouped_column_sums(X, group[a:b], n_groups)
        gene_mean = sums[:len(fit_types)].sum(axis=0) / int(fit.sum())
        mismatch = type_mismatch_genes(sums, fit_types, datasets, src.genes)
        bins, gene_bin, used = make_bins(src.genes, pos, gene_mean, exclude, bin_genes, min_gene_mean,
                                         exclude_pattern=RIBOSOMAL if exclude_ribosomal else None,
                                         genome_build=genome_build, centromeres=centromeres, max_gap_bp=max_gap_bp)
        used.attrs['gene_namespace'] = gene_namespace
        J = len(bins)
        bins['arm'] = arm_labels(bins, centromeres, genome_build)
        Y = open_memmap(out / 'bin_counts.npy', np.float32, (src.n_cells, J))
        library = np.empty(src.n_cells)
        for a, b, X in src.chunks(chunk):
            Y[a:b], library[a:b] = bin_sums(X, gene_bin, J)
        Y.flush()
        n_genes_input = src.n_genes
        # Per fit reference type: counts of each used gene, mixed with a specimen's fitted reference weights to
        # split each bin's expectation across its genes (log_ratio_gene).
        ref_gene = {key: sums[t][used.column.to_numpy()] for t, key in enumerate(fit_types)}
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
    status[discovery & (Lbin < params.min_umi)] = 'insufficient_depth'
    clone = np.full(len(meta), '', dtype=object)
    post_max = np.full(len(meta), np.nan)
    p_normal = np.full(len(meta), np.nan)
    margin = np.full(len(meta), np.nan)
    confirmation_status = np.where(confirmation, 'not_fitted', 'not_requested').astype(object)
    confirmation_records, confirmation_contrasts = [], []
    programmes, depth_gamma = {}, {}
    post_store, tables = {}, {k: [] for k in ['clones', 'segments', 'profiles', 'differences', 'log', 'tree',
                                             'reference_check', 'specimens', 'allele_segments', 'arms']}
    allele_depth = np.zeros(len(meta))
    in_sample_ref = np.zeros(len(meta), bool)
    allele_snps = np.zeros(len(meta), dtype=int)
    for spec in todo:
        rows = np.flatnonzero(discovery & meta.specimen.eq(spec).to_numpy() & (Lbin >= params.min_umi))
        if len(rows) < params.min_cells:
            status[discovery & meta.specimen.eq(spec).to_numpy() & (Lbin >= params.min_umi)] = 'too_few_cells'
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
                params.n_programmes, arm_starts(bins, centromeres, genome_build) if params.programme_window <= 0 else reset,
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
        mu_spec, type_keys, n_in_sample, cand = mu[:, cols], [types[c] for c in cols], 0, None
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
                  genes=(gene_loader(counts, used, rows),
                         dict(counts=[ref_gene[k] for k in type_keys if k in ref_gene], gbin=used.bin.to_numpy(),
                              n_bins=J, in_sample=cand if n_in_sample else None)))
        tables['specimens'][-1]['in_sample_reference_cells'] = n_in_sample
        # Confirmation eligibility must not condition on held-out test-feature counts.
        # score_confirmation abstains when proposal counts cannot support assignment.
        confirm_rows = np.flatnonzero(confirmation & meta.specimen.eq(spec).to_numpy())
        if len(confirm_rows):
            result = confirm_specimen(spec, rows, confirm_rows, meta, Y, res, tables['differences'], out)
            scored = result.pop('scored')
            names = result.pop('candidate_names')
            candidate = scored['candidate']
            eligible = candidate >= 0
            clone[confirm_rows[eligible]] = np.asarray(names, dtype=object)[candidate[eligible]]
            status[confirm_rows[eligible]] = 'confirmation_candidate'
            status[confirm_rows[~eligible]] = 'confirmation_unresolved'
            confirmation_status[confirm_rows] = np.where(eligible, 'proposal_scored', 'unusable_proposal_evidence')
            confirmation_records.append(result)
            confirmation_contrasts.extend(dict(specimen=spec, **row) for row in scored['contrasts'])
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
    cells_out = meta.assign(library_UMI=library, binned_UMI=Lbin,
                            library_count=library, binned_count=Lbin, count_units='input_count_units',
                            allele_depth=allele_depth, allele_snps=allele_snps,
                            in_sample_reference=in_sample_ref,
                            status=status, clone=clone, confirmation_status=confirmation_status,
                            evidence_status=np.where(np.isin(status, ['assigned', 'confirmation_candidate']),
                                                     'candidate_RNA_group', 'unresolved'),
                            calibration_status='unvalidated',
                            clone_posterior=post_max, p_normal=p_normal, assignment_margin=margin,
                            malignant_call=(query & (status == 'assigned')
                                            & ~pd.Series(clone).str.endswith('_c0').to_numpy()
                                            & (p_normal < 0.5)))
    write_table(cells_out, out / 'cells.tsv.gz')
    write_json(dict(specimens=confirmation_records, requested=bool(confirmation.any()),
                    statistical_independence='not_established', specimen_error_control=False,
                    note='Frozen expression-only candidate evaluation; nominal test diagnostics do not validate DNA clones.'),
               out / 'confirmation.json')
    write_table(pd.DataFrame(confirmation_contrasts, columns=['specimen', 'name', 'candidate_a', 'candidate_b',
                'direction', 'n_cells', 'test_counts', 'test_variance', 'testable', 'log_odds_contrast', 'approximate_se',
                'signed_z', 'nominal_one_sided_p', 'multiple_testing']), out / 'confirmation.tsv')
    for name, rows in tables.items():
        if name == 'allele_segments' and allele_data is None:
            continue
        frame = pd.DataFrame(rows) if rows else pd.DataFrame(columns=EMPTY_TABLE_COLUMNS[name])
        write_table(frame, out / f'{name}.tsv.gz' if name == 'profiles' else out / f'{name}.tsv')
    qc = dict(haplodose=__version__, schema_version=2, started=started, wall_seconds=time.monotonic() - t0,
              max_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2,
              n_cells=len(meta), n_genes_input=n_genes_input, n_genes_used=len(used), n_bins=J,
              n_query=int(query.sum()), n_fit_reference=int(fit.sum()), n_holdout_reference=int(hold.sum()),
              n_discovery_query=int(discovery.sum()), n_confirmation_query=int(confirmation.sum()),
              n_confirmation_scored=int((status == 'confirmation_candidate').sum()),
              n_confirmation_unresolved=int((confirmation & (status != 'confirmation_candidate')).sum()),
              reference_groups=types, specimens_fitted=len(post_store),
              programme_prefilter_selected_svd_variance_share={d: None if p is None else p.prefilter_selected_svd_variance_share
                                                              for d, p in programmes.items()},
              alleles=None if allele_data is None else dict(
                  path=str(alleles), n_snps=int(allele_data.A.shape[1]), excluded_rows=allele_data.excluded,
                  exclude_genes=str(allele_exclude_genes), exclude_sites=str(allele_exclude_sites),
                  params=asdict(allele_params)),
              native_kernels=native_available(),
              device=device, backend=backend, devices=[str(d) for d in jax.devices()],
              host=platform.node(), python=sys.version.split()[0], versions=versions(),
              params=dict(asdict(params), bin_genes=bin_genes, min_gene_mean=min_gene_mean, chunk=chunk,
                          query_roles=list(query_roles), exclude_genes=str(exclude_genes),
                          n_excluded_names=len(exclude), exclude_ribosomal=exclude_ribosomal,
                          discovery_fraction=discovery_fraction),
              n_type_mismatch_genes=len(mismatch), type_mismatch_policy='diagnostic_only',
              feature_selection='fit_reference_only', calibration_status='unvalidated',
              coordinate_system='zero-based-half-open', genome_build=genome_build,
              legacy_fields={'malignant_call': 'exploratory reference-difference heuristic; not DNA-validated'},
              note='Relative RNA dosage states and RNA-inferred clones; not absolute DNA copy number '
                   'or validated DNA clones.')
    write_json(qc, out / 'qc.json')
    return qc


def analysis_roles(meta, query, params, discovery_fraction):
    """Assign roles from IDs/metadata before reading query counts or computing query diagnostics."""
    meta = meta.copy()
    explicit = 'analysis_split' in meta
    if explicit and discovery_fraction is not None:
        raise ValueError('Use explicit analysis_split metadata or discovery_fraction, not both')
    if explicit:
        if not meta.loc[query, 'analysis_split'].isin(['discovery', 'confirmation']).all():
            raise ValueError('Query analysis_split must be discovery or confirmation')
    else:
        meta['analysis_split'] = np.where(query, 'discovery', 'not_query')
    for spec in sorted(meta.loc[query, 'specimen'].unique()):
        rows = np.flatnonzero(query & meta.specimen.eq(spec).to_numpy())
        if discovery_fraction is not None:
            roles = assign_roles(meta.cell_id.to_numpy()[rows], discovery_fraction=discovery_fraction,
                                 seed=params.seed, min_discovery=params.min_cells, min_confirmation=1)
            meta.loc[rows, 'analysis_split'] = np.where(roles['discovery'], 'discovery', 'confirmation')
        if meta.analysis_split.iloc[rows].eq('confirmation').any():
            if meta.analysis_split.iloc[rows].eq('discovery').sum() < params.min_cells:
                raise ValueError(f'Specimen {spec} has fewer than min_cells discovery cells')
    meta.loc[~query, 'analysis_split'] = 'not_query'
    return meta


def confirm_specimen(spec, discovery_rows, confirmation_rows, meta, Y, res, differences, out):
    """Freeze discovery arrays and candidates; confirmation never updates the fitted result or tables."""
    model = res['model']
    order = clone_order(res['R'])
    names = [f'{spec}_c{i}' for i in range(len(order))]
    prop = np.arange(Y.shape[1]) % 2 == 0
    test = ~prop
    contrasts = []
    for row in differences:
        if row['specimen'] != spec:
            continue
        region = tuple(range(row['first_bin'], row['last_bin'] + 1))
        mask = np.isin(np.arange(Y.shape[1]), region)
        if not (test & mask).any() or not (test & ~mask).any():
            continue
        contrasts.append(Contrast(names.index(row['clone_a']), names.index(row['clone_b']), region,
                                  1 if row['direction'] == 'higher' else -1, f'contrast_{len(contrasts)}'))
    # Proposal-depth centring is determined solely from discovery cells. It approximates the original
    # full-depth nuisance relation; no confirmation depth or counts feed back into the fitted model.
    totals = np.concatenate([np.asarray(Y[discovery_rows[a:a + 4096]], np.float64)[:, prop].sum(axis=1)
                             for a in range(0, len(discovery_rows), 4096)])
    frozen = FrozenDiscovery(log_rates=model.log_rates(res['doses'])[order], phi=model.phi,
                             proposal_mask=prop, test_mask=test, contrasts=tuple(contrasts),
                             discovery_cell_ids=tuple(meta.cell_id.to_numpy()[discovery_rows]),
                             depth_gamma=model.gamma, proposal_log_depth_centre=float(np.log(np.maximum(totals, 1.)).mean()),
                             programme_loadings=None if model.programmes is None else model.W,
                             programme_precision=None if model.programmes is None else model.programmes.prior_prec,
                             candidate_prior=res['pi'][order])
    key = hashlib.sha256(str(spec).encode()).hexdigest()[:16]
    artifact_path = out / f'confirmation_artifact_{key}.json'
    payload = asdict(frozen)
    payload.update(specimen=spec, candidate_names=names, channel='expression_only',
                   bins_sha256=file_record(out / 'bins.tsv')['sha256'],
                   discovery_assignment_sha256=hashlib.sha256(np.ascontiguousarray(res['R']).tobytes()).hexdigest(),
                   coordinate_system='zero-based-half-open', feature_roles='even_proposal_odd_test',
                   model_scope='proposal-normalized composite NB; discovery fits may use all discovery features')
    write_json(payload, artifact_path)
    candidates, accum = [], [dict(counts=np.zeros((2, 2)), variance=np.zeros((2, 2)), n=np.zeros(2, int))
                            for _ in contrasts]
    for a in range(0, len(confirmation_rows), 4096):
        block = confirmation_rows[a:a + 4096]
        scored = score_confirmation(np.asarray(Y[block], np.float64), meta.cell_id.to_numpy()[block], frozen)
        candidates.append(scored['candidate'])
        for total, row in zip(accum, scored['contrasts']):
            total['counts'] += row['test_counts']
            total['variance'] += row['test_variance']
            total['n'] += row['n_cells']
    scored = dict(candidate=np.concatenate(candidates), evidence=scored['evidence'],
                  contrasts=[contrast_diagnostic(c, t['counts'], t['variance'], t['n'].tolist())
                             for c, t in zip(contrasts, accum)])
    return dict(specimen=spec, n_discovery=len(discovery_rows), n_confirmation=len(confirmation_rows),
                artifact=file_record(artifact_path), evidence=scored['evidence'],
                channel='expression_only', n_frozen_contrasts=len(contrasts), scored=scored, candidate_names=names)


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


def gene_shares(ref_counts, weights, gbin, n_bins, prior=0.5):
    """Each used gene's share of its bin in the specimen's fitted reference (shares sum to 1 per bin).

    ref_counts [T, G]: used-gene counts per reference type, in the order of the model's mixture weights [T].
    Gene rates per type (counts + prior over the type's total) are mixed with those weights, as the bin
    rates are, so a bin's expectation and its gene shares describe the same reference."""
    C = np.atleast_2d(np.asarray(ref_counts, np.float64)) + prior
    rate = np.asarray(weights, np.float64) @ (C / C.sum(axis=1, keepdims=True))
    tot = np.bincount(gbin, weights=rate, minlength=n_bins)
    return rate / tot[gbin]


def gene_level_log_ratio(y, E_bin, share, gbin, r_bin, min_expected=1.0):
    """Gene-level bin log ratio (log_ratio_gene): median over the bin's genes of log((y_g + 0.5) / (E_g + 0.5)).

    E_g splits the bin's model expectation by reference gene shares. Summed-count bins let a few highly
    expressed genes dominate (Kinker single-line rho vs DepMap 0.852 summed vs 0.862 per-gene median,
    results/singleline_decomp_20261101_2354). Genes with E_g < min_expected are left out (their log ratio
    shrinks to 0 in small clones). The result is shifted so that
    its median over the bins with eligible genes equals r_bin's median over those same bins (r_bin is the
    summed-count, centred profile); bins without eligible genes keep r_bin. Reported next to
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
        with CountSource(counts, gene_namespace=used.attrs.get('gene_namespace', 'custom')) as src:
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
        load, mix = genes
        gbin = mix['gbin']
        C = [np.asarray(c, np.float64) for c in mix['counts']]
        if mix['in_sample'] is not None:                   # in-sample normal type: its own cells' gene counts
            Gi, _ = load(np.asarray(mix['in_sample']))
            C.append(np.asarray(Gi.sum(axis=0), np.float64).ravel())
        assert len(C) == len(model.weights), (len(C), len(model.weights))
        share = gene_shares(np.vstack(C), model.weights, gbin, mix['n_bins'])
        Yg = 0.0
        for a in range(0, len(rows), chunk):
            G, _ = load(np.arange(a, min(a + chunk, len(rows))))
            Yg = Yg + np.asarray(G.T @ R[a:a + chunk].astype(np.float64)).T
    chrom = bins.chromosome.to_numpy()
    nb = len(bins)
    al = model.alleles
    arm = bins.arm.to_numpy() if 'arm' in bins else arm_labels(bins)
    for k, name in enumerate(names):
        segs = segments_table(prof['post'][k], prof['path'][k], prof['dose'][k], Yk[k], prof['E'][k], reset,
                              model.beta, prof['tau'][k], model.pp, reference_var=model.v)
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
        # log_dose: this clone's own HMM fit (diagnostic, as segments.tsv); assignment_dose: the dose that set the
        # final cell likelihoods (shared-profile dose from fit_specimen; differs from log_dose for tumour clones).
        assign = np.asarray(res['doses'])[order[k]]
        tables['profiles'].extend(dict(specimen=spec, clone=name, bin=j, log_ratio=float(r[j]),
                                       log_ratio_gene=float(rg[j]),
                                       log_dose=float(prof['dose'][k][j]), assignment_dose=float(assign[j]),
                                       cn_state=int(prof['path'][k][j] + 1))
                                  for j in range(nb))
        overlap = np.nan if k == 0 or model.programmes is None else model.programmes.overlap(prof['dose'][k],
                                                                                             model.base)
        tables['clones'].append(dict(specimen=spec, clone=name, column=k, normal_like=k == 0,
                                     n_cells_map=int((mapc == k).sum()), n_eff=float(R[:, k].sum()),
                                     pooled_UMI=float(Lk[k]), pooled_count=float(Lk[k]),
                                     count_units='input_count_units', resolved_segments=resolved,
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
            allele = (Ar, Dr, model.B, model.alleles.params.cell_shrink, model.alleles.params.cell_kappa)
            Ad = Dr.copy()
            ck = model.alleles.params.cell_kappa
            probability = 0.5 if ck is None else rng.beta(ck / 2, ck / 2, size=len(Dr.data))
            Ad.data = rng.binomial(Dr.data.astype(np.int64), probability).astype(np.float32)
            draw_allele = (Ad, Dr, model.B, model.alleles.params.cell_shrink, ck)
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
