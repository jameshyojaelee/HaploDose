"""Command line: `haplodose run` and `haplodose simulate`."""
from __future__ import annotations

import argparse
import json

from .clones import CloneParams
from .profile import ProfileParams


def main(argv=None):
    ap = argparse.ArgumentParser(prog='haplodose', description=__doc__)
    sub = ap.add_subparsers(dest='command', required=True)
    r = sub.add_parser('run', help='infer relative RNA profiles and candidate CNA groups')
    r.add_argument('--counts', required=True, help='10x H5, h5ad (CSR X) or CSR H5')
    r.add_argument('--cells', required=True, help='TSV in matrix order: cell_id, specimen, patient, role, ref_group, ref_split')
    r.add_argument('--positions', required=True, help='TSV: gene, chromosome, start, end; non-symbol ribosomal filtering requires gene_symbol')
    r.add_argument('--out', required=True, help='new or empty output directory')
    r.add_argument('--exclude-genes', help='one gene name per line, never used as features')
    r.add_argument('--alleles', help='Numbat-style phased allele table (cell, CHROM, POS, cM, AD, DP, GT, gene)')
    r.add_argument('--allele-genome-build', help='required with --alleles; must match the expression annotation build')
    r.add_argument('--allele-exclude-genes', help='genes whose SNPs are never used (e.g. imprinted genes)')
    r.add_argument('--allele-exclude-sites', help='sites never used, one chrN:pos per line')
    r.add_argument('--keep-ribosomal', action='store_true',
                   help='use ribosomal protein genes as features (excluded by default)')
    r.add_argument('--bin-genes', type=int, default=10)
    r.add_argument('--min-gene-mean', type=float, default=0.02)
    r.add_argument('--chunk', type=int, default=4096)
    r.add_argument('--query-roles', default='tumour_epithelium,query')
    r.add_argument('--specimens', nargs='*')
    r.add_argument('--min-cells', type=int, default=30)
    r.add_argument('--event-windows', type=int, nargs='*', default=[],
                   help='sliding-window split-proposal sizes in proposal bins; none given = arms only')
    r.add_argument('--min-umi', type=float, default=100.0)
    r.add_argument('--leiden-resolution', type=float, default=1.0)
    r.add_argument('--p-change', type=float, default=2e-3)
    r.add_argument('--min-diff-bins', type=int, default=3)
    r.add_argument('--programmes', type=int, default=10,
                   help='normal-cell expression programmes used as per-cell nuisance terms (0: off)')
    r.add_argument('--programme-window', type=int, default=0,
                   help='bins of running mean removed from programme loadings (0: arm means, hg38 centromeres)')
    r.add_argument('--depth-source', choices=['query', 'reference'], default='query',
                   help='cells the per-bin depth slopes are estimated from; reference: fit reference cells, '
                        'held fixed (use when depth tracks tumour content, e.g. Visium HD bins)')
    r.add_argument('--seed', type=int, default=20261031)
    r.add_argument('--split-mode', choices=['evenodd', 'thin'], default='evenodd',
                   help='how clone-split proposals and their tests are kept apart: even vs odd bins (default) or '
                        'negative-binomial count splitting of every bin (thin)')
    r.add_argument('--split-null-model', choices=['nb', 'cell_region', 'cov'], default='nb',
                   help='null for --split-null-draws: fitted NB only, plus per-cell, per-chromosome-half rate effects '
                        '(cell_region), or those effects with their full covariance across regions (cov)')
    r.add_argument('--split-null-initial', action='store_true',
                   help='with --split-null-draws: pool the initial tumour groups into one clone (splits rebuild them)')
    r.add_argument('--split-null-draws', type=int, default=0,
                   help='experimental selected-parent simulation rank against this many one-group draws; '
                        'not a calibrated workflow error rate (0: off)')
    r.add_argument('--discovery-fraction', type=float,
                   help='experimental seeded cell split for discovery and independent-cell confirmation; '
                        'alternatively supply analysis_split in metadata')
    r.add_argument('--in-sample-reference', action='store_true',
                   help='experimental opt-in: add allele-balanced query cells to the reference')
    r.add_argument('--device', choices=['cpu', 'gpu', 'auto'], default='cpu',
                   help='JAX kernels on cpu (default), gpu (needs the cuda12 extra and a visible GPU) or auto')
    b = sub.add_parser('bulk', help='relative RNA profile per bulk sample against normal samples')
    b.add_argument('--counts', required=True, help='TSV genes x samples, first column gene name')
    b.add_argument('--samples', required=True, help='TSV: sample, role (tumour | normal)')
    b.add_argument('--positions', required=True, help='TSV: gene, chromosome, start, end; non-symbol ribosomal filtering requires gene_symbol')
    b.add_argument('--out', required=True, help='new or empty output directory')
    b.add_argument('--exclude-genes', help='one gene name per line, never used as features')
    b.add_argument('--bin-genes', type=int, default=10)
    b.add_argument('--device', choices=['cpu', 'gpu', 'auto'], default='cpu',
                   help='JAX kernels on cpu (default), gpu or auto')
    b.add_argument('--min-cpm', type=float, default=1.0, help='minimum mean CPM in normals')
    b.add_argument('--programmes', type=int, default=10, help='normal-sample expression programmes (0: off)')
    b.add_argument('--bin-stat', choices=['sum', 'median'], default='sum',
                   help='bin value: pooled counts, or median of per-gene log ratios to the normals')
    b.add_argument('--rna-profile-smoothing', choices=['none', 'predictive', 'native-regions', 'posterior-regions'], default='none',
                   help='optional descriptive RNA profile: predictive smoothing, native regions or boundary averaging (region modes require sum)')
    for parser in (r, b):
        parser.add_argument('--counts-kind', choices=['raw'], required=True,
                            help='explicit raw-count declaration; normalized/estimated values are unsupported')
        parser.add_argument('--genome-build', required=True,
                            help='annotation assembly; hg38 bundled, others require --centromeres')
        parser.add_argument('--coordinate-system', choices=['one-based-closed', 'zero-based-half-open'],
                            required=True, help='positions input convention (internal output is zero-based half-open)')
        parser.add_argument('--gene-namespace', choices=['symbol', 'ensembl', 'custom'], required=True)
        parser.add_argument('--centromeres', help='TSV chromosome,p_arm_end; zero-based boundary for declared assembly')
        parser.add_argument('--max-gap-bp', type=int,
                            help='experimental: restart bins/HMM across larger genomic gaps (default: off)')
    s = sub.add_parser('simulate', help='write a synthetic dataset with known clones')
    s.add_argument('--out', required=True)
    s.add_argument('--seed', type=int, default=20261031)
    s.add_argument('--depth', type=float, default=3000.0)
    a = ap.parse_args(argv)
    if a.command == 'simulate':
        from .simulate import write_dataset
        print(write_dataset(a.out, seed=a.seed, depth=a.depth))
        return
    if a.command == 'bulk':
        from .bulk import run_bulk
        print(run_bulk(a.counts, a.samples, a.positions, a.out, exclude_genes=a.exclude_genes,
                       bin_genes=a.bin_genes, min_cpm=a.min_cpm, n_programmes=a.programmes,
                       bin_stat=a.bin_stat, device=a.device, counts_kind=a.counts_kind,
                       genome_build=a.genome_build, coordinate_system=a.coordinate_system,
                       gene_namespace=a.gene_namespace, centromeres=a.centromeres, max_gap_bp=a.max_gap_bp,
                       rna_profile_smoothing=a.rna_profile_smoothing))
        return
    from .pipeline import run
    params = CloneParams(min_cells=a.min_cells, min_umi=a.min_umi, leiden_resolution=a.leiden_resolution,
                         n_programmes=a.programmes, event_windows=tuple(a.event_windows),
                         programme_window=a.programme_window, depth_source=a.depth_source, seed=a.seed,
                         split_mode=a.split_mode, split_null_draws=a.split_null_draws,
                         split_null_model=a.split_null_model, split_null_initial=a.split_null_initial,
                         in_sample_reference=a.in_sample_reference,
                         profile=ProfileParams(p_change=a.p_change, min_bins=a.min_diff_bins))
    qc = run(a.counts, a.cells, a.positions, a.out, exclude_genes=a.exclude_genes, bin_genes=a.bin_genes,
             exclude_ribosomal=not a.keep_ribosomal, alleles=a.alleles, allele_genome_build=a.allele_genome_build,
             allele_exclude_genes=a.allele_exclude_genes, allele_exclude_sites=a.allele_exclude_sites,
             min_gene_mean=a.min_gene_mean, chunk=a.chunk, query_roles=tuple(a.query_roles.split(',')),
             specimens=a.specimens, params=params, device=a.device, counts_kind=a.counts_kind,
             genome_build=a.genome_build, coordinate_system=a.coordinate_system,
             gene_namespace=a.gene_namespace, centromeres=a.centromeres, max_gap_bp=a.max_gap_bp,
             discovery_fraction=a.discovery_fraction)
    print(json.dumps({k: qc[k] for k in ['n_cells', 'n_bins', 'specimens_fitted', 'wall_seconds',
                                         'max_rss_gib', 'backend']}))


if __name__ == '__main__':
    main()
