#!/usr/bin/env python3
"""Haplodose input contract for each Andor 2020 gastric line (GSE142750 10x scRNA), read from the GEO tar.

The lines have no normal cells. Reference choice: for line L, the reference is cells from the other
lines, one ref_group per other line (fit/holdout split 50/50). This reference is only valid for
within-line contrasts (differences between L's clones), because the reference expression baseline
of another line differs by its own copy number and transcriptional state. Absolute line profiles
against it are a contrast "L vs the other lines", so the scorer compares them only with the same
contrast in DNA (scDNA consensus of L minus the mean of the other lines' consensus).

Per line: up to --max-query cells of each library of L (role query, specimen = library, patient = L)
and up to --ref-per-line cells per other line. Cells with < --min-umi UMIs are dropped first.
GEO genes files are 2-column (Ensembl, symbol; 33,694 genes) for 7 libraries and 1-column (symbol,
filtered, ~13-14k genes) for 4; genes common to all libraries are used.
Gene positions: GRCh38-2020-A positions.tsv from the CRC YCC230 input, matched by symbol (the GEO
genes files carry symbols only, from an older GENCODE; unmatched genes are dropped by every tool).
Outputs: OUT/<line>/{counts.h5, cells.tsv, positions.tsv}, OUT/manifest.tsv.
"""
import argparse
import gzip
import io
import json
import tarfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import io as sio
from scipy import sparse

TAR = '/gpfs/commons/home/jameslee/Haplodose/data/external/andor2020_gse142750/GSE142750_RAW.tar'
POSITIONS = '/gpfs/commons/home/jameslee/CRC/results/cnv_caller_benchmark/YCC230/input/positions.tsv'
SEED = 20261101


def read_library(tar, prefix):
    get = lambda s: gzip.decompress(tar.extractfile(f'{prefix}_{s}').read())
    genes = pd.read_csv(io.BytesIO(get('genes.tsv.gz')), sep='\t', header=None).iloc[:, -1].astype(str).to_numpy()  # symbol is the last column (1- or 2-column files)
    bcs = pd.read_csv(io.BytesIO(get('barcodes.tsv.gz')), sep='\t', header=None)[0].astype(str).to_numpy()
    X = sparse.csr_matrix(sio.mmread(io.BytesIO(get('matrix.mtx.gz'))).T, dtype=np.float32)  # cells x genes
    return genes, bcs, X


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', required=True)
    ap.add_argument('--max-query', type=int, default=3000)
    ap.add_argument('--ref-per-line', type=int, default=300)
    ap.add_argument('--min-umi', type=float, default=1000)
    args = ap.parse_args()
    rng = np.random.default_rng(SEED)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    tar = tarfile.open(TAR)
    prefixes = sorted({m.name.rsplit('_', 1)[0] for m in tar.getmembers() if m.name.endswith('_matrix.mtx.gz')})
    libs = {}
    for p in prefixes:
        gsm, rest = p.split('_', 1)
        line = rest.rsplit('_', 1)[0] if rest[-2:] in ('_1', '_2') else rest
        genes, bcs, X = read_library(tar, p)
        umi = np.asarray(X.sum(axis=1)).ravel()
        keep = umi >= args.min_umi
        libs[p] = dict(gsm=gsm, line=line, library=rest, genes=genes, cells=np.array([f'{rest}_{b}' for b in bcs])[keep],
                       X=X[keep], n_raw=len(bcs))
    gene_sets = [set(v['genes']) for v in libs.values()]
    common = sorted(set.intersection(*gene_sets))
    pos = pd.read_csv(POSITIONS, sep='\t')
    pos = pos[pos.iloc[:, 0].isin(common)]
    genes = pos.iloc[:, 0].to_numpy()
    for v in libs.values():
        first = ~pd.Index(v['genes']).duplicated()          # duplicated symbols: first occurrence
        idx = np.flatnonzero(first)[pd.Index(v['genes'][first]).get_indexer(genes)]
        v['X'] = v['X'][:, idx]
    lines = sorted({v['line'] for v in libs.values()})
    manifest = []
    for L in lines:
        parts, meta = [], []
        for p, v in libs.items():
            n = v['X'].shape[0]
            if v['line'] == L:
                take = np.sort(rng.choice(n, min(n, args.max_query), replace=False))
                role, spec, group, split = 'query', v['library'], '', ''
                parts.append(v['X'][take])
                meta.append(pd.DataFrame({'cell_id': v['cells'][take], 'specimen': spec, 'patient': L, 'role': role,
                                          'ref_group': group, 'ref_split': split, 'library': v['library']}))
        others = [l for l in lines if l != L]
        for O in others:
            pool = [(p, i) for p, v in libs.items() if v['line'] == O for i in range(v['X'].shape[0])]
            take = rng.choice(len(pool), min(len(pool), args.ref_per_line), replace=False)
            chosen = [pool[i] for i in take]
            split = np.where(rng.random(len(chosen)) < 0.5, 'fit', 'holdout')
            for p in sorted({c[0] for c in chosen}):
                rows = np.array([i for (q, i) in chosen if q == p])
                sp = split[[k for k, c in enumerate(chosen) if c[0] == p]]
                parts.append(libs[p]['X'][rows])
                meta.append(pd.DataFrame({'cell_id': libs[p]['cells'][rows], 'specimen': f'ref_{O}', 'patient': O,
                                          'role': 'reference', 'ref_group': O, 'ref_split': sp,
                                          'library': libs[p]['library']}))
        X = sparse.vstack(parts).tocsr()
        cells = pd.concat(meta, ignore_index=True)
        d = out / L
        d.mkdir()
        with h5py.File(d / 'counts.h5', 'w') as h:
            h.create_dataset('X/data', data=X.data.astype(np.int32), compression='gzip')
            h.create_dataset('X/indices', data=X.indices.astype(np.int32), compression='gzip')
            h.create_dataset('X/indptr', data=X.indptr.astype(np.int64))
            h.create_dataset('genes', data=genes.astype(object), dtype=h5py.string_dtype())
            h.create_dataset('cell_ids', data=cells.cell_id.to_numpy(dtype=object), dtype=h5py.string_dtype())
        cells.to_csv(d / 'cells.tsv', sep='\t', index=False)
        pos.to_csv(d / 'positions.tsv', sep='\t', index=False)
        manifest.append(dict(line=L, n_query=int((cells.role == 'query').sum()), n_reference=int((cells.role == 'reference').sum()),
                             libraries=','.join(sorted(cells.loc[cells.role == 'query', 'library'].unique())), n_genes=len(genes)))
    pd.DataFrame(manifest).to_csv(out / 'manifest.tsv', sep='\t', index=False)
    (out / 'params.json').write_text(json.dumps(dict(vars(args), seed=SEED, tar=TAR, positions=POSITIONS,
                                                     libraries={p: dict(line=v['line'], n_raw=v['n_raw'], n_kept=int(v['X'].shape[0])) for p, v in libs.items()}), indent=2) + '\n')
    print(pd.DataFrame(manifest).to_string())


if __name__ == '__main__':
    main()
