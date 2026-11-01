"""Scale benchmark: stream a large synthetic dataset to a CSR H5, then run Haplodose and record time and memory.

Usage: scale.py OUT [--cells N] [--specimens S] [--simulate-only]

Genes: 22 chromosomes x 500 (about the expressed genes of real data). Each specimen: 10% normal-like cells, clone 1 with two
arm-level events, clone 2 = clone 1 plus one more. Counts are negative binomial at
depth ~1500 UMIs (lognormal). Written in 50k-cell chunks (uint16 data, int16
indices), so memory stays bounded; truth is in cells.tsv (truth_clone).
"""
import argparse
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from haplodose.io import guard, new_output_dir, write_json, write_table

N_CHROM, GPC, DEPTH, PHI, BETA = 22, 500, 1500.0, 0.1, 0.6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--cells', type=int, default=1_000_000)
    ap.add_argument('--specimens', type=int, default=10)
    ap.add_argument('--simulate-only', action='store_true')
    a = ap.parse_args()
    root = new_output_dir(a.out)
    rng = np.random.default_rng(20261101)
    G = N_CHROM * GPC
    chrom = np.repeat([f'chr{c}' for c in range(1, N_CHROM + 1)], GPC)
    within = np.tile(np.arange(GPC), N_CHROM)
    genes = np.array([f'G{i:05d}' for i in range(G)])
    write_table(pd.DataFrame({'gene': genes, 'chromosome': chrom, 'start': within * 1_000_000 + 1,
                              'end': within * 1_000_000 + 50_000}), root / 'positions.tsv')
    base = rng.lognormal(-1.0, 1.2, G)
    types = [base * rng.lognormal(0, 0.6 if t else 0.0, G) for t in range(3)]
    types = [t / t.sum() for t in types]

    populations = []                       # (label, n, rate, meta)
    for t in range(3):
        for split in ('fit', 'holdout'):
            populations.append((f'ref{t}{split}', 2000, types[t],
                                dict(specimen='REF', patient='REF', role='reference', ref_group=f'type{t}',
                                     ref_split=split, truth_clone='reference')))
    per = (a.cells - 12000) // a.specimens
    events = []
    for s in range(a.specimens):
        name = f'S{s:02d}'
        mismatch = types[0] * rng.lognormal(0, 0.15, G)
        ev = rng.choice(N_CHROM, 3, replace=False) + 1
        clone1 = [(f'chr{ev[0]}', 0.0, 1.0, 3), (f'chr{ev[1]}', 0.0, 0.5, 1)]
        clone2 = clone1 + [(f'chr{ev[2]}', 0.5, 1.0, 3)]
        meta = dict(specimen=name, patient=name, role='tumour_epithelium', ref_group='NA', ref_split='NA')
        n_norm = per // 10
        populations.append((f'{name}n', n_norm, types[0], dict(meta, truth_clone='normal')))
        for k, evs in ((1, clone1), (2, clone2)):
            rate = mismatch.copy()
            for c, lo, hi, cn in evs:
                hit = (chrom == c) & (within / GPC >= lo) & (within / GPC < hi)
                rate[hit] *= (cn / 2) ** BETA
                events.append(dict(specimen=name, clone=f'clone{k}', chromosome=c, start_fraction=lo,
                                   end_fraction=hi, cn=cn))
            populations.append((f'{name}c{k}', (per - n_norm) // 2, rate / rate.sum(), dict(meta, truth_clone=f'clone{k}')))
    write_table(pd.DataFrame(events), root / 'events.tsv')

    t0 = time.monotonic()
    meta_rows, n_total = [], 0
    with h5py.File(guard(root / 'counts.h5'), 'w') as h:
        data = h.create_dataset('X/data', shape=(0,), maxshape=(None,), dtype=np.uint16, chunks=(1 << 22,))
        idx = h.create_dataset('X/indices', shape=(0,), maxshape=(None,), dtype=np.int16, chunks=(1 << 22,))
        indptr = [0]
        for label, n, rate, meta in populations:
            for a0 in range(0, n, 50_000):
                m = min(50_000, n - a0)
                L = rng.lognormal(np.log(DEPTH), 0.4, m)
                X = rng.poisson(L[:, None] * rate[None] * rng.gamma(1 / PHI, PHI, (m, G))).astype(np.uint16)
                r, c = np.nonzero(X)
                k = len(data)
                data.resize((k + len(r),))
                idx.resize((k + len(r),))
                data[k:] = X[r, c]
                idx[k:] = c.astype(np.int16)        # < 32768 genes
                indptr.extend((k + np.cumsum(np.bincount(r, minlength=m))).tolist())
                meta_rows.extend(dict(cell_id=f'{label}_{a0 + i}', dataset='scale', **meta) for i in range(m))
                n_total += m
        h.create_dataset('X/indptr', data=np.asarray(indptr, np.int64))
        h.create_dataset('genes', data=genes.astype(object), dtype=h5py.string_dtype())
        h.create_dataset('cell_ids', data=np.array([m['cell_id'] for m in meta_rows], dtype=object),
                         dtype=h5py.string_dtype())
    write_table(pd.DataFrame(meta_rows), root / 'cells.tsv')
    sim = dict(cells=n_total, specimens=a.specimens, simulate_seconds=time.monotonic() - t0,
               simulate_max_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2)
    print(json.dumps(sim), flush=True)
    write_json(sim, root / 'simulation.json')
    if a.simulate_only:
        return
    # Separate process, so its peak memory excludes the simulation.
    subprocess.run([sys.executable, '-s', '-m', 'haplodose.cli', 'run', '--counts', str(root / 'counts.h5'),
                    '--cells', str(root / 'cells.tsv'), '--positions', str(root / 'positions.tsv'),
                    '--out', str(root / 'haplodose')], check=True)
    qc = json.loads((root / 'haplodose' / 'qc.json').read_text())
    cells = pd.read_csv(root / 'haplodose' / 'cells.tsv.gz', sep='\t', usecols=['specimen', 'truth_clone', 'clone',
                                                                               'malignant_call', 'status'])
    from sklearn.metrics import adjusted_rand_score
    q = cells[cells.status.eq('assigned')]
    ari = q.groupby('specimen').apply(lambda g: adjusted_rand_score(g.truth_clone, g.clone), include_groups=False)
    tumour = q.truth_clone.ne('normal')
    summary = dict(sim, wall_seconds=qc['wall_seconds'], max_rss_gib=qc['max_rss_gib'], backend=qc['backend'],
                   n_bins=qc['n_bins'], ari_median=float(ari.median()), ari_min=float(ari.min()),
                   malignant_sensitivity=float(q.malignant_call[tumour].mean()),
                   normal_false_malignant=float(q.malignant_call[~tumour].mean()))
    print(json.dumps(summary), flush=True)
    write_json(summary, root / 'scale_summary.json')
    guard(root / 'pip_freeze.txt').write_text(subprocess.run([sys.executable, '-m', 'pip', 'freeze'],
                                                             capture_output=True, text=True).stdout)


if __name__ == '__main__':
    main()
