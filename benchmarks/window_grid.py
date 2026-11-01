"""Programme window grid on simulated data with known clones.

One specimen per dataset: 10% normal epithelium, clone1 (truncal chr1 gain, chr4
half loss) and clone2 = clone1 plus one subclonal chr6 event of L genes. Two
expression programmes (chr2 50 genes, chr3 150 genes) are switched on in a random
fraction of every population, references included. Each dataset is fitted with
programmes off, with programme windows 25 and 51 bins, and with arm means removed.
Reference: 500 fit cells per type, near SMC (500-1000).

Usage: window_grid.py OUT [--quick]
"""
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from haplodose.clones import CloneParams
from haplodose.io import guard, new_output_dir, write_table
from haplodose.pipeline import run
from haplodose.simulate import write_dataset

GENES_PER_CHROM = 1000
EVENT_START = 0.2
PROGRAMMES = [dict(chromosome='chr2', start=0.40, end=0.45, fold=2.5, frac=0.3),
              dict(chromosome='chr3', start=0.60, end=0.75, fold=1.8, frac=0.25)]
SETTINGS = {'off': CloneParams(n_programmes=0), 'w25': CloneParams(programme_window=25),
            'w51': CloneParams(programme_window=51), 'arm': CloneParams(programme_window=0)}


def coords(chromosome, lo, hi):
    # simulate.py places gene i of a chromosome at i * 100 kb.
    return chromosome, lo * GENES_PER_CHROM * 100_000, hi * GENES_PER_CHROM * 100_000


def overlaps(frame, region):
    c, a, b = region
    return frame.chromosome.eq(c) & (frame.start < b) & (frame.end > a)


def score(out, length):
    cells = pd.read_csv(out / 'cells.tsv.gz', sep='\t')
    cells = cells[cells.status.ne('not_query')]
    clones = pd.read_csv(out / 'clones.tsv', sep='\t')
    tumour = set(clones.loc[~clones.normal_like.astype(bool), 'clone'])
    truth_event = coords('chr6', EVENT_START, EVENT_START + length / GENES_PER_CHROM)
    prog_regions = [coords(p['chromosome'], p['start'], p['end']) for p in PROGRAMMES]
    diff_path = out / 'differences.tsv'
    diff = pd.read_csv(diff_path, sep='\t') if diff_path.stat().st_size > 1 else pd.DataFrame(
        columns=['clone_a', 'clone_b', 'chromosome', 'start', 'end'])
    between = diff[diff.clone_a.isin(tumour) & diff.clone_b.isin(tumour)]
    hit = overlaps(between, truth_event)
    seg = pd.read_csv(out / 'segments.tsv', sep='\t')
    seg = seg[seg.clone.isin(tumour) & seg.resolved.astype(bool) & seg.state.ne('neutral')]
    prog_seg = np.zeros(len(seg), bool)
    prog_diff = np.zeros(len(between), bool)
    for region in prog_regions:
        prog_seg |= overlaps(seg, region).to_numpy()
        prog_diff |= overlaps(between, region).to_numpy()
    malignant = cells[cells.truth_clone.isin(['clone1', 'clone2'])]
    normal = cells[cells.truth_clone.eq('normal')]
    ref = pd.read_csv(out / 'reference_check.tsv', sep='\t').iloc[0]
    return dict(n_tumour_clones=len(tumour),
                ari_malignant=adjusted_rand_score(malignant.truth_clone, malignant.clone.fillna('none')),
                ari_all=adjusted_rand_score(cells.truth_clone, cells.clone.fillna('none')),
                subclone_event_found=bool(hit.any()),
                event_genes_called=int(between.loc[hit, 'n_genes'].sum()),
                false_tumour_differences=int((~hit).sum()),
                programme_differences=int(prog_diff.sum()),
                programme_segments=int(prog_seg.sum()),
                malignant_sensitivity=float(malignant.malignant_call.astype(bool).mean()),
                normal_false_malignant=float(normal.malignant_call.astype(bool).mean()),
                ref_observed=float(ref.fraction_called_aneuploid), ref_expected=float(ref.expected_fraction),
                max_programme_overlap=float(clones.programme_overlap.max()))


def main():
    root = new_output_dir(sys.argv[1])
    quick = '--quick' in sys.argv
    lengths = [60] if quick else [30, 60, 120, 250, 500]
    seeds = [1] if quick else [1, 2]
    directions = {'gain': 3} if quick else {'gain': 3, 'loss': 1}
    rows = []
    for length in lengths:
        for direction, cn in directions.items():
            for seed in seeds:
                name = f'L{length}_{direction}_s{seed}'
                spec = [dict(name='G', n_cells=800, normal_frac=0.1, clones=[
                    dict(frac=0.55, events=[('chr1', 0.0, 1.0, 3), ('chr4', 0.0, 0.5, 1)]),
                    dict(frac=0.45, events=[('chr1', 0.0, 1.0, 3), ('chr4', 0.0, 0.5, 1),
                                            ('chr6', EVENT_START, EVENT_START + length / GENES_PER_CHROM, cn)])])]
                data = write_dataset(root / 'data' / name, seed=20261031 + seed, n_chrom=8, ref_cells_per_type=1000,
                                     genes_per_chrom=GENES_PER_CHROM, specimens=spec, programmes=PROGRAMMES)
                for setting, params in SETTINGS.items():
                    out = root / 'runs' / name / setting
                    t0 = time.time()
                    run(data / 'counts.h5', data / 'cells.tsv', data / 'positions.tsv', out, params=params)
                    row = dict(dataset=name, length=length, direction=direction, seed=seed, setting=setting,
                               wall_s=round(time.time() - t0, 1), **score(out, length))
                    rows.append(row)
                    print(json.dumps(row), flush=True)
                    write_table(pd.DataFrame(rows), root / 'grid.tsv')
    grid = pd.DataFrame(rows)
    summary = grid.groupby(['setting', 'length']).agg(
        n=('dataset', 'size'), event_found=('subclone_event_found', 'mean'), ari=('ari_malignant', 'mean'),
        tumour_clones=('n_tumour_clones', 'mean'), false_diffs=('false_tumour_differences', 'sum'),
        prog_diffs=('programme_differences', 'sum'), prog_segs=('programme_segments', 'sum'),
        ref_obs=('ref_observed', 'mean'), ref_exp=('ref_expected', 'mean')).reset_index()
    write_table(summary, root / 'summary.tsv')
    print(summary.to_string(index=False))
    freeze = subprocess.run([sys.executable, '-m', 'pip', 'freeze'], capture_output=True, text=True).stdout
    guard(root / 'pip_freeze.txt').write_text(freeze)


if __name__ == '__main__':
    main()
