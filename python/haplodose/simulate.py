"""Synthetic scRNA counts with known copy-number clones.

Gene rates follow reference cell types; tumour clones share one epithelial
baseline with a specimen-level gene mismatch (sd tau) and differ by contiguous
dosage events, expression scaled by (cn / 2) ** beta. Counts are gamma-Poisson
(negative binomial). Normal-like epithelial cells use the epithelial reference
baseline without mismatch. Depth changes composition: each gene's share scales
with (L / depth) ** gamma_g, gamma_g tied to expression level and chromosome, so
shallow and deep cells differ chromosome-wide without any copy-number change.
A specimen's `depth_mix` [(fraction, multiplier), ...] makes depth bimodal.
`programmes` [dict(chromosome, start, end, fold, frac), ...] switch on a contiguous
gene cluster (fold change) in a random fraction of every cell population,
reference cells included: expression that mimics copy-number change.

With `alleles=True`, heterozygous SNPs (rate `snp_rate` per gene) get
haplotype-specific copies per event (loss (1,0), gain (2,1), amplification
(3,1), and copy-neutral LOH (2,0) for events written (chrom, lo, hi, 2, 'loh')).
Reads covering a SNP are a binomial share of the gene's counts; the B-haplotype
share follows the copies, with a per-SNP allele-specific expression offset. The
reported phase flips at random switch errors (one per `switch_every` SNPs on
average), so naive pooling of oriented counts cancels real imbalance. Reference
cells take patient `ref_patient` (set it to the specimen's patient for alleles).

With `spatial=True`, specimen cells get x, y: clone k fills the band
k - 1 <= x < k, normal-like cells are spread over all bands.
"""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from .io import guard, new_output_dir, write_table

DEFAULT_SPECIMENS = [
    dict(name='S1', n_cells=600, normal_frac=0.15, clones=[
        dict(frac=0.6, events=[('chr3', 0.0, 1.0, 3), ('chr8', 0.5, 1.0, 1)]),
        dict(frac=0.4, events=[('chr3', 0.0, 1.0, 3), ('chr8', 0.5, 1.0, 1), ('chr12', 0.0, 0.5, 3)]),
    ]),
    dict(name='S2', n_cells=400, normal_frac=0.0, clones=[
        dict(frac=1.0, events=[('chr7', 0.0, 1.0, 3), ('chr18', 0.0, 1.0, 1)]),
    ]),
    dict(name='S3', n_cells=400, normal_frac=0.0, depth_mix=[(0.5, 0.25), (0.5, 1.0)], clones=[
        dict(frac=1.0, events=[('chr5', 0.0, 1.0, 3), ('chr14', 0.0, 1.0, 1)]),
    ]),
]


def simulate(seed=20261031, n_chrom=22, genes_per_chrom=200, n_types=3, ref_cells_per_type=200,
             specimens=DEFAULT_SPECIMENS, depth=3000.0, phi=0.1, beta=0.6, tau=0.15, type_sd=0.6,
             depth_slope=0.15, extra_genes=('IGF2',), programmes=(), alleles=False, snp_rate=0.5,
             snp_capture=0.05, switch_every=20, ase_sd=0.2, ref_patient='REF', spatial=False):
    rng = np.random.default_rng(seed)
    chrom = np.repeat([f'chr{c}' for c in range(1, n_chrom + 1)], genes_per_chrom)
    within = np.tile(np.arange(genes_per_chrom), n_chrom)
    genes = [f'G{i:05d}' for i in range(len(chrom))]
    # Extra named genes (e.g. a masked gene) sit at the end of chr11 at high expression.
    for name in extra_genes:
        genes.append(name)
        chrom = np.append(chrom, 'chr11')
        within = np.append(within, genes_per_chrom)
    G = len(genes)
    positions = pd.DataFrame({'gene': genes, 'chromosome': chrom, 'start': within * 100_000 + 1,
                              'end': within * 100_000 + 50_000})
    base = rng.lognormal(-1.0, 1.2, G)
    base[len(genes) - len(extra_genes):] = base.max()
    types = [base * rng.lognormal(0, type_sd if t else 0.0, G) for t in range(n_types)]
    types = [t / t.sum() for t in types]       # type 0 is the epithelial reference baseline
    frac_pos = within / genes_per_chrom
    z = np.log(base)
    chrom_effect = dict(zip(np.unique(chrom), rng.normal(0, 0.1, len(np.unique(chrom)))))
    gamma = depth_slope * ((z - z.mean()) / z.std() + rng.normal(0, 0.3, G)) + np.vectorize(chrom_effect.get)(chrom)
    rows, meta, allele_rows = [], [], []
    if alleles:
        snp_gene = np.repeat(np.arange(G), rng.poisson(snp_rate, G))
        snp_gene = snp_gene[snp_gene < G - len(extra_genes)]
        snp_pos = (positions.start.to_numpy()[snp_gene] + rng.integers(0, 50_000, len(snp_gene)))
        order = np.lexsort((snp_pos, np.vectorize(lambda c: int(c[3:]))(chrom[snp_gene])))
        snp_gene, snp_pos = snp_gene[order], snp_pos[order]
        true_b_alt = rng.random(len(snp_gene)) < 0.5            # ALT allele sits on haplotype B
        flips = rng.random(len(snp_gene)) < 1.0 / switch_every
        new_chrom = np.r_[True, chrom[snp_gene][1:] != chrom[snp_gene][:-1]]
        block = np.zeros(len(snp_gene), bool)
        for i in range(1, len(snp_gene)):
            block[i] = False if new_chrom[i] else block[i - 1] ^ flips[i]
        reported_b_alt = true_b_alt ^ block
        ase = rng.normal(0, ase_sd, len(snp_gene))

    def draw_alleles(X, copies_b, label_rows):
        """copies_b [G] B-haplotype share per gene for this population."""
        d = rng.binomial(X[:, snp_gene].astype(np.int64), snp_capture)
        b = copies_b[snp_gene]
        inner = (b > 0) & (b < 1)
        logit = np.log(np.where(inner, b, 0.5) / (1 - np.where(inner, b, 0.5))) + ase
        b = np.where(inner, 1 / (1 + np.exp(-logit)), b)
        x = rng.binomial(d, b[None])
        r, c = np.nonzero(d)
        alt = np.where(true_b_alt[c], x[r, c], d[r, c] - x[r, c])
        allele_rows.append(pd.DataFrame(dict(cell=np.asarray(label_rows)[r], snp_id=c,
                                             CHROM=[ch[3:] for ch in chrom[snp_gene][c]], POS=snp_pos[c],
                                             cM=snp_pos[c] / 1e6, REF='A', ALT='C', AD=alt, DP=d[r, c],
                                             GT=np.where(reported_b_alt[c], '0|1', '1|0'),
                                             gene=np.asarray(genes)[snp_gene][c])))

    def draw(rate, n, label, depth_mix=((1.0, 1.0),), copies_b=None, xy=None, **cols):
        mult = np.concatenate([np.full(int(round(n * f)), m) for f, m in depth_mix])[:n]
        mult = np.concatenate([mult, np.full(n - len(mult), depth_mix[-1][1])])
        L = rng.lognormal(np.log(depth), 0.4, n) * mult
        share = rate[None] * (L[:, None] / depth) ** gamma[None]
        for prog in programmes:
            on = rng.random(n) < prog['frac']
            hit = (chrom == prog['chromosome']) & (frac_pos >= prog['start']) & (frac_pos < prog['end'])
            share[np.ix_(on, hit)] *= prog['fold']
        share /= share.sum(axis=1, keepdims=True)
        lam = L[:, None] * share * rng.gamma(1 / phi, phi, (n, G))
        X = rng.poisson(lam).astype(np.int32)
        rows.append(sparse.csr_matrix(X))
        ids = [f'{label}_{len(meta) + i}' for i in range(n)]
        for i in range(n):
            meta.append(dict(cell_id=ids[i], **cols) if xy is None else dict(cell_id=ids[i], x=xy[i, 0], y=xy[i, 1], **cols))
        if alleles:
            draw_alleles(X, np.full(G, 0.5) if copies_b is None else copies_b, ids)

    for t in range(n_types):
        for split in ('fit', 'holdout'):
            draw(types[t], ref_cells_per_type // 2, f'ref{t}', specimen='REF', patient=ref_patient, dataset='sim',
                 role='reference', ref_group=f'type{t}', ref_split=split, truth_clone='reference')
    events_out = []
    for spec in specimens:
        mismatch = types[0] * rng.lognormal(0, tau, G)
        n_norm = int(round(spec['n_cells'] * spec['normal_frac']))
        n_bands = len(spec['clones'])
        if n_norm:
            xy = np.column_stack([rng.random(n_norm) * n_bands, rng.random(n_norm)]) if spatial else None
            draw(types[0], n_norm, spec['name'], xy=xy, specimen=spec['name'], patient=spec['name'], dataset='sim',
                 role='tumour_epithelium', ref_group='NA', ref_split='NA', truth_clone='normal')
        remaining = spec['n_cells'] - n_norm
        for k, cl in enumerate(spec['clones'], start=1):
            n = int(round(remaining * cl['frac']))
            rate = mismatch.copy()
            copies_b = np.full(G, 0.5)
            for event in cl['events']:
                c, lo, hi, cn = event[:4]
                loh = len(event) > 4 and event[4] == 'loh'
                hit = (chrom == c) & (frac_pos >= lo) & (frac_pos < hi)
                rate[hit] *= (cn / 2.0) ** beta
                major = {1: 1, 2: 2 if loh else 1, 3: 2, 4: 3}[cn]
                # Which haplotype carries the major copies is fixed per event (seeded by its position).
                b_major = np.random.default_rng([seed, int(c[3:]), int(lo * 1000)]).random() < 0.5
                share_major = major / cn
                copies_b[hit] = share_major if b_major else 1 - share_major
                events_out.append(dict(specimen=spec['name'], clone=f'clone{k}', chromosome=c,
                                       start_fraction=lo, end_fraction=hi, cn=cn, loh=loh, n_genes=int(hit.sum())))
            xy = np.column_stack([k - 1 + rng.random(n), rng.random(n)]) if spatial else None
            draw(rate / rate.sum(), n, spec['name'], depth_mix=spec.get('depth_mix', ((1.0, 1.0),)), copies_b=copies_b, xy=xy,
                 specimen=spec['name'], patient=spec['name'], dataset='sim', role='tumour_epithelium',
                 ref_group='NA', ref_split='NA', truth_clone=f'clone{k}')
    X = sparse.vstack(rows).tocsr()
    allele_table = pd.concat(allele_rows, ignore_index=True) if alleles else None
    return X, genes, pd.DataFrame(meta), positions, pd.DataFrame(events_out), allele_table


def write_dataset(out, **kwargs):
    out = new_output_dir(out)
    X, genes, cells, positions, events, allele_table = simulate(**kwargs)
    with h5py.File(guard(out / 'counts.h5'), 'w') as h:
        h.create_dataset('X/data', data=X.data.astype(np.int32), compression='gzip')
        h.create_dataset('X/indices', data=X.indices.astype(np.int32), compression='gzip')
        h.create_dataset('X/indptr', data=X.indptr.astype(np.int64))
        h.create_dataset('genes', data=np.asarray(genes, dtype=object), dtype=h5py.string_dtype())
        h.create_dataset('cell_ids', data=cells.cell_id.to_numpy(dtype=object), dtype=h5py.string_dtype())
    write_table(cells, out / 'cells.tsv')
    write_table(positions, out / 'positions.tsv')
    write_table(events, out / 'events.tsv')
    if allele_table is not None:
        write_table(allele_table, out / 'alleles.tsv.gz')
    return Path(out)
