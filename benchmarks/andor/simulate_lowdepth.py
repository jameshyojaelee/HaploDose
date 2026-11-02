#!/usr/bin/env python3
"""Simulate the SRA cell-tagged low-depth subset for a line from its full-depth per-cell bins.

Why: the SNU-601 SRA run kept no cell read groups (all reads carry the 8-mer sample index), so the
low-depth pipeline cannot be run on real SNU-601 SRA data. To measure what the low-depth truth loses,
each full-depth cell is binomially thinned to a fragment count drawn from the cell-tagged per-cell
fragment counts of the other lines, with the per-bin cell-tagged bias of those lines applied
(median over lines of cell-tagged share / index share per 500 kb bin; shared by cells of a line).
The full-depth pooled counts are written as one index read group ("BULK") for the bulk profile.

Usage: simulate_lowdepth.py FULL_BINS_DIR SRA_BINS_DIR LINE OUT_BINS_DIR
Writes OUT_BINS_DIR/LINE.tsv.gz, LINE.rg.txt (copied), sim_params.json. Seeded.
"""
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20261101
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]


def read_bins(path):
    df = pd.read_csv(path, sep='\t', header=None, names=['rg', 'chromosome', 'bin', 'n'], dtype={'rg': str})
    return df[df.chromosome.isin(AUTOSOMES)]


def main(full_dir, sra_dir, line, out):
    full_dir, sra_dir, out = Path(full_dir), Path(sra_dir), Path(out)
    out.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(SEED)

    # Per-cell fragment counts and per-bin bias from the other lines' SRA cell-tagged subsets.
    frags, ratios, used = [], [], []
    for p in sorted(sra_dir.glob('*.tsv.gz')):
        other = p.name[:-len('.tsv.gz')]
        if other == line:
            continue
        df = read_bins(p)
        cell = df.rg.str.contains('-')
        if cell.sum() == 0:
            continue
        used.append(other)
        frags.append(df[cell].groupby('rg').n.sum().to_numpy())
        c = df[cell].groupby(['chromosome', 'bin']).n.sum()
        i = df[~cell].groupby(['chromosome', 'bin']).n.sum()
        j = pd.concat([c.rename('c'), i.rename('i')], axis=1).fillna(0)
        j = j[(j.c > 0) & (j.i > 0)]
        ratios.append(np.log2((j.c / j.c.sum()) / (j.i / j.i.sum())).rename(other))
    frags = np.concatenate(frags)
    bias = 2 ** pd.concat(ratios, axis=1).median(axis=1)

    full = read_bins(full_dir / f'{line}.tsv.gz')
    w = bias.reindex(pd.MultiIndex.from_frame(full[['chromosome', 'bin']])).fillna(1.0).to_numpy()
    full = full.assign(w=w)
    tot_w = full.assign(nw=full.n * full.w).groupby('rg').nw.sum()
    cells = tot_w.index.to_numpy()
    target = pd.Series(rng.choice(frags, size=len(cells), replace=True), index=cells)
    p = np.minimum(target.reindex(full.rg).to_numpy() * full.w.to_numpy() / tot_w.reindex(full.rg).to_numpy(), 1.0)
    thin = full.assign(n=rng.binomial(full.n.to_numpy(), p))
    thin = thin[thin.n > 0][['rg', 'chromosome', 'bin', 'n']]
    bulk = full.groupby(['chromosome', 'bin']).n.sum().reset_index().assign(rg='BULK')[['rg', 'chromosome', 'bin', 'n']]
    pd.concat([thin, bulk]).to_csv(out / f'{line}.tsv.gz', sep='\t', header=False, index=False)
    shutil.copy(full_dir / f'{line}.rg.txt', out / f'{line}.rg.txt')
    got = thin.groupby('rg').n.sum()
    params = dict(seed=SEED, line=line, source_lines=used, n_cells=int(len(cells)),
                  donor_cells=int(len(frags)), donor_median_fragments=float(np.median(frags)),
                  full_median_fragments=float(full.groupby('rg').n.sum().median()),
                  sim_median_fragments=float(got.median()),
                  bias_bins=int(len(bias)), bias_log2_sd=float(np.log2(bias).std()))
    (out / 'sim_params.json').write_text(json.dumps(params, indent=1))
    print(json.dumps(params, indent=1))


if __name__ == '__main__':
    main(*sys.argv[1:])
