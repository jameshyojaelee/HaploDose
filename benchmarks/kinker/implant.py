"""Kinker subclone design: implant a known copy-number change into a random subset of one line's cells.

A single-line specimen is one population, so it cannot test clone separation. Here, in each rotation of the
v2 design (design.py), every single-line query with >= MIN_CELLS cells gets one event in a random fraction
of its cells. The two groups share the line's expression state, so only the event separates them.
- Setting per line from a 12-cell grid (size x direction x fraction), balanced over all queries (seeded):
  size = 50, 150 or 500 expressed genes (all-cell mean >= 0.05 UMI) contiguous on one autosome;
  loss = in-region counts thinned binomially to 0.5 (2 -> 1 copies);
  gain (GAIN_MODE 'depth_neutral', design v3, default) = every cell of the specimen thinned to 2/3, except the
  region in implanted cells (2 -> 3 copies, relative); both groups keep the same library size, as real gains
  roughly do. GAIN_MODE 'thin_outside' (design v2, dirs up to 20261101_1854) thinned only implanted cells
  outside the region, so they had a third fewer UMIs; tools that model depth as a covariate (Haplodose depth
  slopes and depth-matched pooling) then read the gain as depth.
  fraction = 0.15 or 0.4 of the line's usable cells.
- Only query cells of that rotation are changed, so each rotation has its own counts file.
Outputs: counts_r{0,1}.h5, cells_r{0,1}.tsv (+ implanted column), implants.tsv, positions.tsv link.
Usage: implant.py INPUT_DIR NEW_DIR [thin_outside|depth_neutral]
"""
import itertools
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

SEED = 20261102
MIN_CELLS = 100
MIN_MEAN = 0.05
SIZES = (50, 150, 500)
DIRECTIONS = ('loss', 'gain')
FRACTIONS = (0.15, 0.4)
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]


def read_h5(path):
    with h5py.File(path) as h:
        X = sparse.csr_matrix((h['X/data'][:], h['X/indices'][:], h['X/indptr'][:]))
        return X, h['genes'][:].astype(str), h['cell_ids'][:].astype(str)


def write_h5(path, X, genes, cells):
    with h5py.File(path, 'w') as h:
        h.create_dataset('genes', data=genes.astype(h5py.string_dtype()))
        h.create_dataset('cell_ids', data=cells.astype(h5py.string_dtype()))
        g = h.create_group('X')
        for k, v in [('data', X.data.astype(np.int32)), ('indices', X.indices.astype(np.int32)), ('indptr', X.indptr)]:
            g.create_dataset(k, data=v, compression='gzip')


def expressed_order(X, genes, inp):
    """Expressed autosomal genes ordered by position: DataFrame gene, col, chromosome, start."""
    pos = pd.read_csv(inp / 'positions.tsv', sep='\t').drop_duplicates('gene')
    pos['chromosome'] = 'chr' + pos.chromosome.astype(str).str.removeprefix('chr')
    mean = np.asarray(X.mean(axis=0)).ravel()
    g = pd.DataFrame({'gene': genes, 'col': np.arange(len(genes)), 'mean': mean})
    g = g.merge(pos[['gene', 'chromosome', 'start']], on='gene')
    g = g[g.chromosome.isin(AUTOSOMES) & (g['mean'] >= MIN_MEAN)]
    return g.sort_values(['chromosome', 'start']).reset_index(drop=True)


def pick_region(rng, expr, size, all_pos):
    chroms = expr.groupby('chromosome').size()
    chrom = rng.choice(chroms.index[chroms >= size + 20])
    e = expr[expr.chromosome.eq(chrom)].reset_index(drop=True)
    i = int(rng.integers(0, len(e) - size + 1))
    start, end = int(e.start.iloc[i]), int(e.start.iloc[i + size - 1])
    cols = all_pos.col[all_pos.chromosome.eq(chrom) & all_pos.start.between(start, end)].to_numpy()
    return chrom, start, end, cols


def main(inp, out, gain_mode='depth_neutral'):
    assert gain_mode in ('thin_outside', 'depth_neutral')
    inp, out = Path(inp).resolve(), Path(out)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'positions.tsv').symlink_to(inp / 'positions.tsv')
    rng = np.random.default_rng(SEED)
    rng_neutral = np.random.default_rng(SEED + 1)   # own stream: v3 keeps v2's regions and implanted cells
    X0, genes, cell_ids = read_h5(inp / 'counts.h5')
    expr = expressed_order(X0, genes, inp)
    pos = pd.read_csv(inp / 'positions.tsv', sep='\t').drop_duplicates('gene')
    pos['chromosome'] = 'chr' + pos.chromosome.astype(str).str.removeprefix('chr')
    all_pos = pd.DataFrame({'gene': genes, 'col': np.arange(len(genes))}).merge(pos[['gene', 'chromosome', 'start']], on='gene')
    row = pd.Series(np.arange(len(cell_ids)), index=cell_ids)
    cells = {r: pd.read_csv(inp / f'cells_r{r}.tsv', sep='\t', dtype={'pool': str}) for r in (0, 1)}
    queries = pd.concat([c[c.design.eq('single_line')].groupby('specimen').size().rename('n').reset_index().assign(rotation=r)
                         for r, c in cells.items()])
    queries = queries[queries.n >= MIN_CELLS].sample(frac=1, random_state=SEED).reset_index(drop=True)
    grid = list(itertools.product(SIZES, DIRECTIONS, FRACTIONS))
    rows = []
    for r, c in cells.items():
        X = X0.copy()
        c['implanted'] = False
        for i, q in queries[queries.rotation.eq(r)].iterrows():
            size, direction, fraction = grid[i % len(grid)]
            chrom, start, end, cols = pick_region(rng, expr, size, all_pos)
            ids = c.cell_id[c.specimen.eq(q.specimen) & c.design.eq('single_line')].to_numpy()
            pick = rng.choice(ids, int(round(fraction * len(ids))), replace=False)
            in_region = np.zeros(X.shape[1], bool)
            in_region[cols] = True
            for k in row[pick].to_numpy():
                a, b = X.indptr[k], X.indptr[k + 1]
                hit = in_region[X.indices[a:b]] if direction == 'loss' else ~in_region[X.indices[a:b]]
                p = 0.5 if direction == 'loss' else 2 / 3
                X.data[a:b][hit] = rng.binomial(X.data[a:b][hit], p)
            if direction == 'gain' and gain_mode == 'depth_neutral':
                for k in row[np.setdiff1d(ids, pick)].to_numpy():
                    a, b = X.indptr[k], X.indptr[k + 1]
                    X.data[a:b] = rng_neutral.binomial(X.data[a:b], 2 / 3)
            c.loc[c.cell_id.isin(pick), 'implanted'] = True
            rows.append(dict(rotation=r, specimen=q.specimen, n_cells=len(ids), n_implanted=len(pick),
                             size=size, direction=direction, fraction=fraction, gain_mode=gain_mode, chromosome=chrom,
                             start=start, end=end, n_genes_in_region=len(cols)))
        X.eliminate_zeros()
        write_h5(out / f'counts_r{r}.h5', X, genes, cell_ids)
        c.to_csv(out / f'cells_r{r}.tsv', sep='\t', index=False)
        print(f'rotation {r}: {int(c.implanted.sum())} implanted cells', flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(out / 'implants.tsv', sep='\t', index=False)
    print(t.groupby(['size', 'direction', 'fraction']).size().to_string())


if __name__ == '__main__':
    main(*sys.argv[1:])
