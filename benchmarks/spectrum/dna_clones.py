#!/usr/bin/env python3
"""SPECTRUM DLP+ clone copy-number tables and the DNA part of the claim 4 patient rule (docs/STATUS.md, claim 4).

Input (Synapse syn66366960 download, data/external/spectrum_ov_scdna_synapse):
  <patient>_cell_info.csv.gz   cell_id, leaf_id (the published tumour clone of each tumour cell), n_wgd, tumor_site
  signals_<patient>.h5         AnnData: layers/state = HMMcopy integer total copy number, cells x 500 kb bins;
                               var chr/start/end are GRCh37 (chr1 ends at 249.5 Mb, chromosomes '1'..'22','X')
Clone CN per bin = median layers/state over the clone's cells (as benchmarks/andor/scdna_cn.py clone_bins.tsv).
leaf_id rule: leaf_id is a clone label, read as text. A compound label such as '0/2' or '0/1/3' is ONE published
clone: <patient>_cluster_info.csv.gz lists it as a single leaf with its own cluster_size, and
<patient>_branch_info.csv.gz gives it its own terminal branch (clone_0/2) next to the internal branches, so each cell
has exactly one clone (it is not a cell assigned to several leaves, and not an internal node). It is kept as is,
with '/' written as '_' in column names (clone_0_2). Cells with an empty leaf_id have no published clone and are
dropped. Only cells in cell_info are used; normal cells are not in cell_info.

Patient rule, DNA part (fixed before any RNA tool runs): >= 2 DLP+ tumour clones and >= 1 arm where two clones
differ on >= 30% of its bins by >= 1 copy (score.py dna_truth, hg19). Clones are taken as published, with no size
filter (the rule counts published clones). A patient whose truth flags every arm has no defined AUROC (score.py) and
is not eligible (reason 'all arms differ'; e.g. a WGD and a non-WGD clone differ by >= 1 total copy everywhere). The RNA part of the rule (an unsorted or CD45-negative tumour scRNA
sample) is applied by the RNA sample manifest, not here.

Usage: dna_clones.py SYNAPSE_DIR OUT
Writes OUT/<patient>/clone_bins.tsv (chromosome, start, end, clone_<leaf_id>...), OUT/<patient>/dna_clone_differences.tsv,
OUT/clones.tsv (patient, clone, n_cells, n_wgd) and OUT/patients.tsv (patient, dna_clones, dna_subclonal_arms,
arms, eligible_dna, reason).
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score as sp  # noqa: E402


def decode(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a])


def categorical(f, key):
    categories, codes = decode(f[f'{key}/categories'][:]), f[f'{key}/codes'][:]
    if np.any((codes < 0) | (codes >= len(categories))):
        raise ValueError(f'missing or invalid categorical code: {key}')
    return categories[codes]


def clone_bins(h5, cell_info):
    with h5py.File(h5, 'r') as f:
        var = pd.DataFrame({'chromosome': categorical(f, 'var/chr'), 'start': f['var/start'][:], 'end': f['var/end'][:]})
        cid = pd.Index(decode(f['obs/cell_id'][:]))
        ci = pd.read_csv(cell_info, dtype={'leaf_id': str, 'cell_id': str})
        if cid.hasnans or not cid.is_unique or ci.cell_id.isna().any() or ci.cell_id.duplicated().any():
            raise ValueError('DNA cell identifiers must be unique and nonmissing')
        # leaf_id is a label ('0', '2/4' = published merged leaves); cells with none have no published clone, not used
        ci = ci[ci.leaf_id.notna()].copy()
        labels = ci.leaf_id.drop_duplicates()
        if labels.str.replace('/', '_', regex=False).duplicated().any():
            raise ValueError('DNA clone label sanitization would merge distinct published labels')
        ci['leaf_id'] = ci.leaf_id.str.replace('/', '_', regex=False)
        if ci.empty:
            raise ValueError('no DNA cells have a published clone label')
        idx = cid.get_indexer(ci.cell_id)
        if (idx < 0).any():
            raise ValueError(f'{(idx < 0).sum()} cell_info cells missing from {h5}')
        order = np.argsort(idx)
        state = f['layers/state'][np.sort(idx), :]            # h5py fancy indexing needs increasing indices
    leaf = ci.leaf_id.to_numpy()[order]
    out = var.copy()
    sizes = {}
    for k in sorted(pd.unique(leaf)):
        out[f'clone_{k}'] = np.median(state[leaf == k], axis=0)
        sizes[f'clone_{k}'] = int((leaf == k).sum())
    wgd = ci.groupby('leaf_id').n_wgd.agg(lambda x: x.mode().iloc[0])
    return out, sizes, {f'clone_{k}': int(v) for k, v in wgd.items()}


def main(src, out):
    src, out = Path(src), Path(out)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(out)
    out.mkdir(parents=True, exist_ok=True)
    rows, crow = [], []
    for ci in sorted(src.glob('*_cell_info.csv.gz')):
        p = ci.name.replace('_cell_info.csv.gz', '')
        h5 = src / f'signals_{p}.h5'
        d = out / p
        d.mkdir()
        if not h5.exists():
            rows.append(dict(patient=p, eligible_dna=False, reason='no signals h5'))
            continue
        cb, sizes, wgd = clone_bins(h5, ci)
        cb.to_csv(d / 'clone_bins.tsv', sep='\t', index=False)
        names, all_arms, diff_arms, diffs, qc = sp.dna_truth(d / 'clone_bins.tsv', 'hg19', d)
        diffs.to_csv(d / 'dna_clone_differences.tsv', sep='\t', index=False)
        crow += [dict(patient=p, clone=c, n_cells=sizes[c], n_wgd=wgd[c]) for c in names]
        ok = len(names) >= 2 and 1 <= len(diff_arms) < len(all_arms)
        reason = 'ok' if ok else ('< 2 clones' if len(names) < 2 else 'no differing arm' if not diff_arms else 'all arms differ')
        rows.append(dict(patient=p, dna_clones=len(names), dna_cells=sum(sizes.values()), min_clone_cells=min(sizes.values()),
                         dna_subclonal_arms=len(diff_arms), arms=len(all_arms), eligible_dna=ok, reason=reason, **qc))
        print(rows[-1], flush=True)
    pd.DataFrame(crow).to_csv(out / 'clones.tsv', sep='\t', index=False)
    t = pd.DataFrame(rows)
    t.to_csv(out / 'patients.tsv', sep='\t', index=False)
    print(t.to_string())
    print(f'eligible (DNA part): {int(t.eligible_dna.sum())}/{len(t)}')


if __name__ == '__main__':
    main(*sys.argv[1:])
