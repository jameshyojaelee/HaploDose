#!/usr/bin/env python3
"""Frozen T3 prototype comparator (~/CRC/src/43_cnv_caller/caller.py, read only).

Run with the prototype's env: ~/micromamba/envs/rapids_singlecell/bin/python -B.
The prototype's writers refuse paths outside ~/CRC, so this wrapper imports caller.py
unchanged and replaces only its write helpers (guard, output_dir, atomic_table,
atomic_text) with ones that write here. Bytecode writing is off so nothing lands in
~/CRC/src. Expression only (no --baf). Query roles other than 'tumour_epithelium' are
renamed to it in a staged cells.tsv, because the prototype only calls that role.

Output: native/ holds the prototype's own files; the rest is the Haplodose layout
(cmp_common.write_layout). The prototype calls malignancy only for tumour_epithelium cells, so
held-out reference cells are never malignant by construction (reference_check is 0).
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True

import argparse  # noqa: E402
import os  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cmp_common as C  # noqa: E402

PROTO = Path('/gpfs/commons/home/jameslee/CRC/src/43_cnv_caller')


def _output_dir(path):
    p = C.guard(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def _atomic_table(frame, path):
    p = C.guard(path)
    _output_dir(p.parent)
    fd, name = tempfile.mkstemp(prefix=f'.{p.name}.', dir=p.parent)
    os.close(fd)
    frame.to_csv(name, sep='\t', index=False, na_rep='NA', compression='gzip' if p.suffix == '.gz' else None)
    os.replace(name, p)


def _atomic_text(path, text):
    p = C.guard(path)
    _output_dir(p.parent)
    fd, name = tempfile.mkstemp(prefix=f'.{p.name}.', dir=p.parent)
    with os.fdopen(fd, 'w') as fh:
        fh.write(text)
    os.replace(name, p)


def load_caller():
    sys.path.insert(0, str(PROTO.parents[0] / '25_multimodal_atlas'))
    sys.path.insert(0, str(PROTO))
    import caller
    caller.guard, caller.output_dir = C.guard, _output_dir
    caller.atomic_table, caller.atomic_text = _atomic_table, _atomic_text
    return caller


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--counts', required=True)
    ap.add_argument('--cells', required=True)
    ap.add_argument('--positions', required=True)
    ap.add_argument('--out-parent', required=True, help='a new <out-parent>/prototype_<timestamp> is created')
    ap.add_argument('--query-roles', default='tumour_epithelium,query')
    ap.add_argument('--bin-genes', type=int, default=50)
    ap.add_argument('--chunk', type=int, default=512)
    args = ap.parse_args(argv)

    t0 = time.monotonic()
    out = C.new_run_dir(args.out_parent, 'prototype')
    C.provenance(out, 'prototype', args, prototype=str(PROTO / 'caller.py'))
    C.pip_freeze(out)
    caller = load_caller()

    cells = C.read_cells(args.cells)
    query_roles = args.query_roles.split(',')
    staged = Path(args.cells)
    if set(query_roles) - {'tumour_epithelium'} and cells.role.isin(query_roles).any():
        full = pd.read_csv(args.cells, sep='\t', dtype=str, keep_default_na=False)
        full.loc[full.role.isin(query_roles), 'role'] = 'tumour_epithelium'
        staged = out / 'staged_cells.tsv'
        C.write_table(full, staged)
    native = out / 'native'
    caller.run(argparse.Namespace(counts=args.counts, cells=str(staged), positions=args.positions,
                                  out=str(native), chunk=args.chunk, bin_genes=args.bin_genes,
                                  spatial_k=0, baf=None))

    res = pd.read_csv(native / 'cells.tsv.gz', sep='\t', dtype=str, keep_default_na=False)
    if not res.cell_id.equals(cells.cell_id):
        raise ValueError('Prototype output order differs from cells.tsv')
    malignant = res.malignant_call.eq('True').to_numpy()
    query = C.query_mask(cells, query_roles)
    hold = (cells.role.eq('reference') & cells.ref_split.eq('holdout')).to_numpy()
    analysed = query | hold
    clone = np.where(malignant, res.clone, np.where(analysed, cells.specimen + '_normal', 'NA')).astype(object)
    clone[analysed & (res.clone == 'unresolved').to_numpy() & ~malignant] = 'unresolved'
    # Haplodose layout: units are prototype bins; log_ratio = clone mean of the prototype's cell
    # profiles minus the held-out reference median (its own centring); cn_state from the clone's
    # HMM prototype for tumour clones (cnLOH counts as neutral dosage), calibrated rule otherwise.
    bins = pd.read_csv(native / 'bins.tsv', sep='\t')
    genes_df = (bins.assign(gene=bins.genes.str.split(';')).explode('gene')[['gene', 'chromosome', 'bin']])
    P = np.load(native / 'profiles.npy', mmap_mode='r')
    lib = pd.to_numeric(res.library_UMI).to_numpy()
    centre = np.median(np.asarray(P[hold & (lib >= 300)]), axis=0)
    hold_ok = hold & (lib >= 300)
    Z = np.asarray(P, dtype=np.float32) - centre
    profiles = C.clone_profiles(Z, cells, clone, hold_ok)
    try:
        native_clones = pd.read_csv(native / 'clones.tsv', sep='\t')
    except pd.errors.EmptyDataError:
        native_clones = pd.DataFrame(columns=['specimen', 'clone', 'states'])
    to_cn = np.array([1, 2, 3, 2])
    for _, r in native_clones.iterrows():
        st = np.asarray(str(r.states).split(','), dtype=int)
        sel = profiles.clone.eq(r.clone).to_numpy()
        if sel.any():
            profiles.loc[sel, 'cn_state'] = np.where(st >= 0, to_cn[np.clip(st, 0, 3)], 2)
    normal_like = [c for c in pd.unique(clone) if str(c).endswith('_normal') or c == 'unresolved']
    C.write_layout(out, cells, clone, malignant, bins[['bin', 'chromosome', 'start', 'end', 'n_genes']],
                   genes_df, profiles, normal_like, hold)
    C.write_json(dict(tool='prototype', n_cells=len(cells), n_malignant=int(malignant.sum()),
                      native_qc=str(native / 'qc.json'),
                      holdout_malignant_rate=float(malignant[hold].mean()) if hold.any() else None,
                      wall_seconds=time.monotonic() - t0), out / 'qc.json')
    print(out)


if __name__ == '__main__':
    main()
