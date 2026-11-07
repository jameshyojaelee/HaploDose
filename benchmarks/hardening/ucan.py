"""Locked U-CAN development mechanisms on verified featureCounts observations.

Run with the installed candidate wheel on a compute node. Source code may be
frozen separately from the live --data-root. The two tumour IDs and all available
normal samples are fixed independently of DNA values. No former evaluation
cohort, model setting, reference donor or smoothing width is selected here.

Example (the width is an explicit, previously developed diagnostic setting):
  python ucan.py --data-root /path/to/Haplodose --out NEWROOT --smooth-bins 5

The runner verifies built counts against the original gene-ID featureCounts
matrix and the locked tumours' per-sample exon tables before any inference.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import gzip
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.ndimage import median_filter
from scipy.stats import pearsonr, spearmanr

LOCKED_SAMPLES = ('CRC.SW.U0001.T', 'CRC.SW.U0008.T')
CRC = Path('/gpfs/commons/home/jameslee/CRC')
RAW_SOURCE = CRC / 'data/interim/full_transcriptome_counts/exon_counts.tsv.gz'
ANNOTATION = CRC / 'data/external/ensembl112_fulltx_features/genes.tsv'
QUANTIFIER = CRC / 'src/03_requant/09_quantify_fulltx_sample.sh'
MASK = CRC / 'src/44_haplodose/masks.txt'
VCF_DIR = CRC / 'data/external/eva_cnv/vcf'


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def file_record(path):
    path = Path(path).resolve()
    return dict(path=str(path), bytes=path.stat().st_size, sha256=sha256(path))


def unique(values, description):
    index = pd.Index(values)
    if index.has_duplicates or index.isna().any() or any(not str(x).strip() for x in index):
        raise ValueError(f'{description} must be nonempty unique identifiers')
    return index


def verify_inputs(inp, raw_source, annotation):
    """Reconstruct the builder's symbol mapping, checking every selected observation."""
    meta = pd.read_csv(inp / 'samples.tsv', sep='\t', dtype=str)
    unique(meta['sample'], 'sample IDs')
    if not meta.role.isin(['normal', 'tumour']).all():
        raise ValueError('unexpected sample role in established U-CAN input')
    if not set(LOCKED_SAMPLES) <= set(meta['sample']):
        raise ValueError('locked U-CAN tumour is absent')
    selected = meta[meta.role.eq('normal') | meta['sample'].isin(LOCKED_SAMPLES)].copy().reset_index(drop=True)
    if selected.loc[selected.role.eq('normal'), 'patient'].duplicated().any():
        raise ValueError('expected one normal RNA sample per biological donor')
    if selected.role.eq('normal').sum() < 3:
        raise ValueError('established normal panel is incomplete')
    names = selected['sample'].to_numpy()
    X = pd.read_csv(inp / 'counts.tsv', sep='\t', index_col=0, usecols=['gene', *names])
    X = X.loc[:, names]
    unique(X.index, 'built gene symbols')
    values = X.to_numpy(float)
    if not np.isfinite(values).all() or np.any(values < 0) or np.any(values != np.floor(values)):
        raise ValueError('U-CAN observations are not finite nonnegative integer counts')
    if np.any(values.sum(axis=0) <= 0):
        raise ValueError('U-CAN input contains an empty sample')
    features = pd.read_csv(annotation, sep='\t', dtype={'contig': str, 'gene_id': str, 'gene_name': str})
    features = features[features.biotype.eq('protein_coding') & features.contig.isin([str(i) for i in range(1, 23)])]
    features = features[~features.gene_name.duplicated(keep=False)].set_index('gene_name')
    if not set(X.index) <= set(features.index):
        raise ValueError('built symbols cannot be reconstructed from the established annotation')
    mapping = features.loc[X.index]
    unique(mapping.gene_id, 'source gene IDs')
    original = pd.read_csv(raw_source, sep='\t', index_col=0, usecols=['gene_id', *names])
    unique(original.index, 'original featureCounts gene IDs')
    if not set(mapping.gene_id) <= set(original.index):
        raise ValueError('a built gene is absent from the original featureCounts matrix')
    expected = original.loc[mapping.gene_id, names].to_numpy(float)
    if not np.array_equal(values, expected):
        raise ValueError('built U-CAN counts differ from original gene-ID counts')
    positions = pd.read_csv(inp / 'positions.tsv', sep='\t', dtype={'chromosome': str})
    unique(positions.gene, 'position gene symbols')
    p = positions.set_index('gene').loc[X.index]
    if not (p.chromosome.to_numpy() == mapping.contig.to_numpy()).all():
        raise ValueError('position chromosomes differ from the source annotation')
    if not np.array_equal(p[['start', 'end']].to_numpy(), mapping[['start', 'end']].to_numpy()):
        raise ValueError('position coordinates differ from Ensembl 112 one-based closed annotation')
    per_sample = []
    for sample in LOCKED_SAMPLES:
        path = raw_source.parent / 'per_sample' / f'{sample}.exon.tsv.gz'
        counts = pd.read_csv(path, sep='\t', index_col=0, usecols=['gene_id', sample])
        unique(counts.index, f'per-sample gene IDs for {sample}')
        if not np.array_equal(counts.loc[mapping.gene_id, sample].to_numpy(float), X[sample].to_numpy(float)):
            raise ValueError(f'aggregate matrix disagrees with per-sample exon output for {sample}')
        per_sample.append(path)
    return X, selected, positions, dict(
        status='verified_raw_featureCounts_observations',
        n_genes=len(X), n_normal_donors=int(selected.role.eq('normal').sum()),
        n_locked_tumours=len(LOCKED_SAMPLES), n_values_verified=int(X.size),
        mapping='unique protein-coding autosomal symbols, no library scaling or count transformation',
        coordinate_system='Ensembl 112 one-based-closed input',
        source_quantifier='featureCounts exon/gene_id; unstranded paired fragments, proper pairs, MAPQ>=20, primary, no duplicates',
        source_versions='legacy featureCounts/samtools runtime versions were not recorded by the source script',
        per_sample_files=[str(p) for p in per_sample])


def canonicalize_observations(counts, positions):
    """After raw-source verification, use the installed wheel's exact symbol identity contract.

    Counts are neither summed nor changed. Ambiguous case-folded identities fail
    through validate_feature_ids instead of being silently collapsed.
    """
    from haplodose.io import validate_feature_ids
    source = counts.index.to_numpy(dtype=str)
    canonical = validate_feature_ids(source, 'symbol')
    result = counts.copy()
    result.index = canonical
    result.index.name = 'gene'
    position = positions.copy()
    position['gene'] = validate_feature_ids(position.gene.to_numpy(dtype=str), 'symbol')
    mapping = pd.DataFrame(dict(source_gene=source, canonical_gene=canonical))
    return result, position, mapping


def read_dna(patient):
    """VCF POS/END become zero-based half-open intervals, matching current package bins."""
    paths = list(VCF_DIR.glob(f'*CRC-SW-{patient}-T*.vcf.gz'))
    if len(paths) != 1:
        raise ValueError(f'expected exactly one existing FACETS VCF for {patient}')
    header, rows = {}, []
    with gzip.open(paths[0], 'rt') as stream:
        for line in stream:
            if line.startswith(('##purity=', '##ploidy=')):
                key, value = line[2:].strip().split('=', 1)
                header[key] = float(value)
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            info = dict(item.split('=', 1) for item in fields[7].split(';') if '=' in item)
            cn = info.get('TCN_EM', '.')
            rows.append(dict(chromosome='chr' + fields[0].removeprefix('chr'), start=int(fields[1]) - 1,
                             end=int(info['END']), tcn=float(cn) if cn not in ('.', 'NA') else np.nan))
    if not np.isfinite(header.get('ploidy', np.nan)) or header['ploidy'] <= 0:
        raise ValueError(f'invalid FACETS ploidy for {patient}')
    seg = pd.DataFrame(rows).sort_values(['chromosome', 'start'])
    for chrom, g in seg.groupby('chromosome'):
        if (g.start < 0).any() or (g.end <= g.start).any() or np.any(g.start.to_numpy()[1:] < g.end.to_numpy()[:-1]):
            raise ValueError(f'invalid or overlapping FACETS intervals for {patient} {chrom}')
    if (seg.tcn.dropna() < 0).any():
        raise ValueError('negative FACETS copy number')
    return seg, header, paths[0]


def dna_at_genes(used, seg, ploidy):
    midpoint = ((used.start.to_numpy() + used.end.to_numpy() - 1) // 2)
    result = np.full(len(used), np.nan)
    for chrom, s in seg.groupby('chromosome', sort=False):
        rows = np.flatnonzero(used.chromosome.eq(chrom).to_numpy())
        if not len(rows):
            continue
        start, end = s.start.to_numpy(), s.end.to_numpy()
        idx = np.searchsorted(start, midpoint[rows], side='right') - 1
        ok = (idx >= 0) & (midpoint[rows] < end[np.maximum(idx, 0)])
        cn = s.tcn.to_numpy()[idx[ok]]
        result[rows[ok]] = np.log(np.maximum(cn, .25) / ploidy)
    return result


def reduce_genes(values, gene_bin, n_bins, require_complete=False):
    result = np.full(n_bins, np.nan)
    for j in range(n_bins):
        vector = values[gene_bin == j]
        if len(vector) and (not require_complete or np.isfinite(vector).all()):
            result[j] = np.nanmedian(vector)
    return result


def smooth(values, blocks, width):
    result = np.full_like(values, np.nan, dtype=float)
    for block in pd.unique(blocks):
        take = blocks == block
        result[take] = median_filter(values[take], size=width, mode='nearest')
    return result


def metrics(x, y):
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[valid], np.asarray(y)[valid]
    if len(x) < 4 or np.ptp(y) == 0:
        return dict(n=len(x), spearman=np.nan, pearson=np.nan)
    return dict(n=len(x), spearman=float(spearmanr(x, y).statistic) if np.ptp(x) else 0.0,
                pearson=float(pearsonr(x, y).statistic) if np.ptp(x) else 0.0)


def installed_package_record(data_root):
    import haplodose
    package = Path(haplodose.__file__).resolve().parent
    if package == data_root / 'python/haplodose' or (data_root / 'python') in package.parents:
        raise ValueError('runner requires the installed candidate wheel, not an editable repository package')
    dist = importlib.metadata.distribution('haplodose')
    direct = dist.read_text('direct_url.json')
    if direct and json.loads(direct).get('dir_info', {}).get('editable'):
        raise ValueError('an editable HaploDose distribution is not an immutable candidate wheel')
    paths = sorted(set(package.rglob('*.py')) | set(package.rglob('*.so')))
    return dict(version=dist.version, package_root=str(package), direct_url=direct,
                files=[file_record(path) for path in paths])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--smooth-bins', type=int, required=True,
                        help='explicit fixed diagnostic width; identical operator for every method')
    parser.add_argument('--device', choices=['cpu', 'gpu'], default='cpu')
    args = parser.parse_args(argv)
    if args.smooth_bins < 1 or args.smooth_bins % 2 == 0:
        raise ValueError('smoothing width must be a positive odd number of current genomic bins')
    data_root = args.data_root.resolve()
    inp = data_root / 'results/ucan_bulk_input_20261101_1352'
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    package = installed_package_record(data_root)
    X, meta, positions, verified = verify_inputs(inp, RAW_SOURCE, ANNOTATION)
    X, positions, identity = canonicalize_observations(X, positions)
    identity.to_csv(out / 'gene_identity_map.tsv', sep='\t', index=False)
    verified['canonical_identity'] = dict(namespace='symbol', operation='installed-wheel strip/uppercase',
                                          n_changed=int(identity.source_gene.ne(identity.canonical_gene).sum()),
                                          ambiguity_policy='reject duplicates; never aggregate observations')
    input_dir = out / 'input'
    input_dir.mkdir()
    X.to_csv(input_dir / 'counts.tsv', sep='\t')
    meta.to_csv(input_dir / 'samples.tsv', sep='\t', index=False)
    positions.to_csv(input_dir / 'positions.tsv', sep='\t', index=False)
    sources = [inp / 'counts.tsv', inp / 'samples.tsv', inp / 'positions.tsv', RAW_SOURCE,
               ANNOTATION, QUANTIFIER, MASK, data_root / 'benchmarks/ucan_bulk/make_input.py',
               *[Path(p) for p in verified['per_sample_files']]]
    provenance = dict(analysis='locked_UCAN_raw_count_mechanisms', evidence_role='development',
                      independent_evaluation=False, biological_n=2, locked_samples=LOCKED_SAMPLES,
                      selection='two locked tumour IDs plus all original normal samples; no DNA selection',
                      source_validation=verified, source_files=[file_record(p) for p in sources],
                      executed_script=file_record(__file__), installed_candidate=package,
                      argv=sys.argv if argv is None else argv, slurm_job=os.environ.get('SLURM_JOB_ID'),
                      hostname=os.uname().nodename, units='natural_log_relative_rate',
                      coordinate_contract='hg38, symbol, one-based-closed input; internal zero-based-half-open',
                      smoothing=dict(width_bins=args.smooth_bins, operator='running_median_nearest',
                                     boundary='current genomic block; never crosses arm or declared gap'),
                      model_parameters=dict(bin_genes=10, min_cpm=1.0, n_programmes=10, bin_stat='sum',
                                            max_gap_bp=None, device=args.device, exclude_genes=str(MASK)))
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2))
    print(json.dumps(dict(source_validation=verified, action='current installed wheel inference starting')), flush=True)
    from haplodose.bulk import BULK_PROFILE, fit_sample, reference_model, run_bulk
    from haplodose.bins import arm_labels, arm_starts, chromosome_starts
    from haplodose.profile import log_ratio
    provenance['model_parameters']['profile'] = asdict(BULK_PROFILE)
    run_bulk(input_dir / 'counts.tsv', input_dir / 'samples.tsv', input_dir / 'positions.tsv',
             out / 'haplodose', exclude_genes=MASK, bin_genes=10, min_cpm=1.0,
             n_programmes=10, params=BULK_PROFILE, bin_stat='sum', device=args.device,
             counts_kind='raw', genome_build='hg38', coordinate_system='one-based-closed',
             gene_namespace='symbol', centromeres=None, max_gap_bp=None)
    hd = out / 'haplodose'
    bins = pd.read_csv(hd / 'bins.tsv', sep='\t')
    used = pd.read_csv(hd / 'genes.tsv.gz', sep='\t')
    samples = pd.read_csv(hd / 'samples.tsv', sep='\t')
    if not np.array_equal(samples['sample'], meta['sample']):
        raise ValueError('candidate bulk output changed sample order')
    G = X.loc[used.gene, meta['sample']].to_numpy(float)
    gene_bin = used.bin.to_numpy(int)
    onehot = sparse.csr_matrix((np.ones(len(used)), (np.arange(len(used)), gene_bin)),
                              shape=(len(used), len(bins)))
    Y = np.asarray((onehot.T @ G).T)
    normal = meta.role.eq('normal').to_numpy()
    reference = reference_model(Y[normal], 10, arm_starts(bins))
    reset = chromosome_starts(bins)
    arms = arm_labels(bins)
    blocks = bins.genomic_block.to_numpy() if 'genomic_block' in bins else bins.chromosome.to_numpy()
    current_lr = np.load(hd / 'log_ratio.npy')
    segments = pd.read_csv(hd / 'segments.tsv', sep='\t')
    unique(used.gene, 'candidate used genes')
    if not np.array_equal(used.start.to_numpy(), positions.set_index('gene').start.reindex(used.gene).to_numpy() - 1):
        raise ValueError('candidate gene coordinates do not obey the declared one-based input conversion')
    score_rows, traces, bin_rows, dna_sources = [], [], [], []
    for sample in LOCKED_SAMPLES:
        i = int(np.flatnonzero(meta['sample'].eq(sample).to_numpy())[0])
        f = fit_sample(Y[i], *reference, reset, BULK_PROFILE)
        pooled_expectation = Y[i].sum() * reference[0]
        pooled_lr = log_ratio(Y[i], pooled_expectation)
        pooled_lr -= np.median(pooled_lr)
        corrected_lr = log_ratio(Y[i], f['E'])
        corrected_lr -= np.median(corrected_lr)
        replay_error = float(np.max(np.abs(corrected_lr - current_lr[i])))
        if not np.allclose(corrected_lr, current_lr[i], atol=1e-5, rtol=1e-6):
            raise ValueError('stage replay does not reproduce the candidate exported log ratio')
        gene_mu = G[:, normal].sum(axis=1) / G[:, normal].sum()
        gene_expectation = G[:, i].sum() * gene_mu
        gene_ratio = log_ratio(G[:, i], gene_expectation)
        gene_ratio -= np.median(gene_ratio)
        # Established naive normalization uses the full input protein-coding
        # autosomal library, then exactly the candidate's feature/bin support.
        full_totals = X.sum(axis=0).to_numpy(float)
        cpm = G / full_totals[None] * 1e6
        logged = np.log(cpm + 1)
        naive = logged[:, i] - logged[:, normal].mean(axis=1)
        naive -= np.median(naive)
        exported = np.full(len(bins), np.nan)
        for row in segments[segments['sample'].eq(sample)].itertuples():
            exported[row.first_bin:row.last_bin + 1] = row.rna_log_effect
        methods = dict(gene_ratio=reduce_genes(gene_ratio, gene_bin, len(bins)),
                       naive=reduce_genes(naive, gene_bin, len(bins)), pooled_bin_ratio=pooled_lr,
                       hd_ratio=corrected_lr, hd_posterior=f['dose'], hd_segment_export=exported)
        for method, values in list(methods.items()):
            methods[method + '_smoothed'] = smooth(values, blocks, args.smooth_bins)
        # DNA is first opened only after reference fitting and all method outputs.
        patient = meta.iloc[i]['patient']
        dna, header, vcf = read_dna(patient)
        dna_sources.append(file_record(vcf))
        truth_genes = dna_at_genes(used, dna, header['ploidy'])
        truth = reduce_genes(truth_genes, gene_bin, len(bins), require_complete=True)
        common = np.isfinite(truth) & np.logical_and.reduce([np.isfinite(v) for v in methods.values()])
        frame = bins.copy()
        frame['sample'], frame['dna_ln_tcn_over_ploidy'], frame['common_callable'] = sample, truth, common
        frame['pooled_counts'], frame['pooled_expected'], frame['programme_adjusted_expected'] = Y[i], pooled_expectation, f['E']
        frame['programme_shift'], frame['reference_variance'] = f['shift'], reference[1]
        for name, values in methods.items():
            frame[name] = values
            score_rows.append(dict(sample=sample, method=name, resolution='common_bin', **metrics(values[common], truth[common])))
            arm_table = pd.DataFrame(dict(arm=arms[common], value=values[common], truth=truth[common]))
            for aggregation in ('median', 'mean'):
                a = arm_table.groupby('arm')[['value', 'truth']].agg(aggregation)
                score_rows.append(dict(sample=sample, method=name, resolution='common_arm_' + aggregation,
                                       **metrics(a.value.to_numpy(), a.truth.to_numpy())))
        frame['dna_gene_coverage'] = [float(np.isfinite(truth_genes[gene_bin == j]).mean()) for j in range(len(bins))]
        bin_rows.append(frame)
        stages = [('gene_to_count_weighted_bin', 'gene_ratio', 'pooled_bin_ratio'),
                  ('normal_programme_adjustment', 'pooled_bin_ratio', 'hd_ratio'),
                  ('HMM_posterior_projection', 'hd_ratio', 'hd_posterior'),
                  ('segment_export_averaging', 'hd_posterior', 'hd_segment_export')]
        stages += [('matched_running_median', name, name + '_smoothed')
                   for name in ('gene_ratio', 'naive', 'pooled_bin_ratio', 'hd_ratio', 'hd_posterior', 'hd_segment_export')]
        for label, before, after in stages:
            delta = methods[after] - methods[before]
            before_score, after_score = metrics(methods[before][common], truth[common]), metrics(methods[after][common], truth[common])
            traces.append(dict(sample=sample, stage=label, before=before, after=after,
                               rms_change=float(np.sqrt(np.mean(delta[common] ** 2))) if common.any() else np.nan,
                               delta_spearman=after_score['spearman'] - before_score['spearman'],
                               delta_pearson=after_score['pearson'] - before_score['pearson'],
                               replay_max_abs_error=replay_error, beta=f['beta'], tau=f['tau'],
                               centre=f['centre'], purity=header.get('purity'), ploidy=header['ploidy']))
        print(f'locked mechanism complete: {sample}, {int(common.sum())}/{len(bins)} common bins', flush=True)
    pd.DataFrame(score_rows).to_csv(out / 'scores.tsv', sep='\t', index=False)
    pd.DataFrame(traces).to_csv(out / 'stage_effects.tsv', sep='\t', index=False)
    pd.concat(bin_rows, ignore_index=True).to_csv(out / 'bin_stage_trace.tsv.gz', sep='\t', index=False)
    provenance['dna_sources'] = dna_sources
    provenance['derived_input_files'] = [file_record(input_dir / name) for name in ('counts.tsv', 'samples.tsv', 'positions.tsv')]
    provenance['candidate_run'] = file_record(hd / 'run.json')
    provenance['elapsed_seconds'] = time.time() - started
    provenance['interpretation'] = ('two selected development examples, no inferential p-values or independent validation; '
                                    'DNA target is natural-log TCN/ploidy at common gene-anchored bins, then identical arm aggregation; '
                                    'neither posterior support nor descriptive arm calls establish calibrated error control')
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2))
    (out / 'summary.json').write_text(json.dumps(dict(
        evidence_role='development', independent_evaluation=False, biological_n=2,
        n_reference_donors=int(normal.sum()), locked_samples=LOCKED_SAMPLES,
        model_parameters=provenance['model_parameters'], scores=score_rows, stages=traces,
        elapsed_seconds=provenance['elapsed_seconds']), indent=2))


if __name__ == '__main__':
    main()
