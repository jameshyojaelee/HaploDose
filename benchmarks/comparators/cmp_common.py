"""Shared input and output helpers for the comparator runners.

Input contract (same as `haplodose run`): counts (10x H5, h5ad with CSR X, or CSR H5
with X/*, genes, cell_ids), cells.tsv (cell_id, specimen, patient, role, ref_group,
ref_split) and positions.tsv (gene, chromosome, start, end).

Common output, one directory per run:
  cells.tsv  cell_id, specimen, clone, malignant_call (every input cell; clone is NA
             for cells the tool did not analyse, e.g. fit reference cells)
  arms.tsv   specimen, clone, arm, n_units, then median_log_ratio and/or state

Only h5py, numpy, pandas and scipy are imported here, so every tool env can use it.
"""
from __future__ import annotations

import datetime as dt
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

REQUIRED_CELL_COLUMNS = ['cell_id', 'specimen', 'patient', 'role', 'ref_group', 'ref_split']
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]
CENTROMERES = Path(__file__).resolve().parents[2] / 'python/haplodose/data/centromeres_hg38.tsv'
# Outputs never go under the CRC project (its source, data and results stay read only here).
FORBIDDEN = [Path('/gpfs/commons/home/jameslee/CRC').resolve()]
SEED = 20261101


def guard(path) -> Path:
    path = Path(os.path.realpath(os.path.abspath(path)))
    for root in FORBIDDEN:
        if path == root or root in path.parents:
            raise PermissionError(f'Refusing to write under {root}: {path}')
    return path


def new_run_dir(parent, tool) -> Path:
    """New timestamped directory <parent>/<tool>_<YYYYmmdd_HHMMSS>; never reused."""
    stamp = dt.datetime.now().strftime('%Y%m%d_%H%M%S')
    out = guard(Path(parent) / f'{tool}_{stamp}')
    out.mkdir(parents=True, exist_ok=False)
    return out


def write_table(frame: pd.DataFrame, path) -> Path:
    path = guard(path)
    tmp = guard(path.with_name('.' + path.name + '.tmp'))
    frame.to_csv(tmp, sep='\t', index=False, na_rep='NA')
    os.replace(tmp, path)
    return path


def write_json(obj, path) -> Path:
    path = guard(path)
    path.write_text(json.dumps(obj, indent=2, default=str) + '\n')
    return path


def pip_freeze(out) -> None:
    res = subprocess.run([sys.executable, '-m', 'pip', 'freeze'], capture_output=True, text=True)
    guard(Path(out) / 'pip_freeze.txt').write_text(res.stdout)


def provenance(out, tool, args, **extra) -> None:
    write_json(dict(tool=tool, args=vars(args) if hasattr(args, '__dict__') else args, seed=SEED,
                    started=dt.datetime.now().isoformat(timespec='seconds'), host=platform.node(),
                    python=sys.executable, slurm_job=os.environ.get('SLURM_JOB_ID'), **extra),
               Path(out) / 'run.json')


def read_cells(path) -> pd.DataFrame:
    cells = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    missing = set(REQUIRED_CELL_COLUMNS) - set(cells.columns)
    if missing:
        raise ValueError(f'cells.tsv lacks columns: {sorted(missing)}')
    if cells.cell_id.duplicated().any():
        raise ValueError('Duplicated cell IDs in cells.tsv')
    return cells


def read_positions(path) -> pd.DataFrame:
    pos = pd.read_csv(path, sep='\t')
    if list(pos.columns[:4]) != ['gene', 'chromosome', 'start', 'end']:
        pos = pos.rename(columns={pos.columns[0]: 'gene'})
    pos['chromosome'] = pos.chromosome.astype(str)
    pos.loc[~pos.chromosome.str.startswith('chr'), 'chromosome'] = 'chr' + pos.chromosome
    return pos.drop_duplicates('gene')


def _h5ad_index(group):
    key = group.attrs.get('_index', '_index')
    return group[key].asstr()[:]


def read_counts(path, cell_ids=None):
    """Cells x genes CSR (float32), gene names, cell ids; reordered to `cell_ids` if given."""
    with h5py.File(path, 'r') as h:
        if 'matrix' in h:
            g = h['matrix']
            genes, cells = g['features/name'].asstr()[:], g['barcodes'].asstr()[:]
        elif 'obs' in h and 'var' in h:
            g = h['X']
            if g.attrs.get('encoding-type') != 'csr_matrix':
                raise ValueError('h5ad X must be stored as csr_matrix')
            genes, cells = _h5ad_index(h['var']), _h5ad_index(h['obs'])
        else:
            g = h['X']
            genes, cells = h['genes'].asstr()[:], h['cell_ids'].asstr()[:]
        X = sparse.csr_matrix((g['data'][:].astype(np.float32), g['indices'][:], g['indptr'][:]),
                              shape=(len(cells), len(genes)))
    if cell_ids is not None:
        idx = pd.Index(cells).get_indexer(pd.Index(cell_ids))
        if (idx < 0).any():
            raise ValueError(f'{int((idx < 0).sum())} cells in cells.tsv are absent from the counts')
        X, cells = X[idx], np.asarray(cell_ids)
    return X, np.asarray(genes), np.asarray(cells)


def collapse_genes(X, genes):
    """Sum columns that share a gene name (10x features can repeat symbols)."""
    uniq, inv = np.unique(genes, return_inverse=True)
    if len(uniq) == len(genes):
        return X, genes
    M = sparse.csr_matrix((np.ones(len(genes), np.float32), (np.arange(len(genes)), inv)),
                          shape=(len(genes), len(uniq)))
    return (X @ M).tocsr(), uniq


def arm_of(chromosome, start) -> np.ndarray:
    """'chrNp'/'chrNq' by start vs the hg38 p-arm end; same rule as haplodose.bins.arm_labels."""
    cen = pd.read_csv(CENTROMERES, sep='\t', comment='#').set_index('chromosome').p_arm_end
    chromosome = pd.Series(np.asarray(chromosome, dtype=str))
    q = np.asarray(start, dtype=float) >= chromosome.map(cen).fillna(np.inf).to_numpy()
    return (chromosome + np.where(q, 'q', 'p')).to_numpy()


def query_mask(cells, query_roles) -> np.ndarray:
    return cells.role.isin(query_roles).to_numpy()


def common_cells(cells, clone, malignant) -> pd.DataFrame:
    return pd.DataFrame({'cell_id': cells.cell_id.to_numpy(), 'specimen': cells.specimen.to_numpy(),
                         'clone': clone, 'malignant_call': np.asarray(malignant, dtype=bool)})


# ---------------------------------------------------------------------------------------------
# Haplodose-shaped output, so benchmarks/ucan_spatial/dna_truth.py and
# benchmarks/cell_lines/evaluate.py score every tool with the same code.
# A "unit" is the tool's genomic unit (infercnvpy window, prototype bin, 50-gene bin).
# log_dose is set equal to log_ratio: comparators have no separate dose estimate.
# cn_state follows Haplodose: 1 loss, 2 neutral, 3 gain.

def gene_bins(genes, pos, size=50):
    """Consecutive `size`-gene bins per autosome, ordered by start. Returns bins, gene table."""
    p = pos.set_index('gene').reindex(pd.Index(genes))
    p = p[p.chromosome.isin(AUTOSOMES)].copy()
    p['order'] = p.chromosome.str[3:].astype(int)
    p = p.sort_values(['order', 'start'])
    p['bin'] = -1
    rows, b = [], 0
    for chrom, g in p.groupby('chromosome', sort=False):
        for lo in range(0, len(g), size):
            block = g.iloc[lo:lo + size]
            if len(block) < max(5, size // 3) and rows and rows[-1]['chromosome'] == chrom:
                p.loc[block.index, 'bin'] = b - 1          # short tail joins the previous bin
                rows[-1]['end'] = int(block.end.max())
                rows[-1]['n_genes'] += len(block)
                continue
            p.loc[block.index, 'bin'] = b
            rows.append(dict(bin=b, chromosome=chrom, start=int(block.start.min()), end=int(block.end.max()),
                             n_genes=len(block)))
            b += 1
    genes_df = p.reset_index().rename(columns={'index': 'gene'})[['gene', 'chromosome', 'start', 'end', 'bin']]
    return pd.DataFrame(rows), genes_df


def calibrated_states(mean, n, hold, z=3.0, q=0.99):
    """cn_state per unit for a clone mean profile from n cells, calibrated on held-out reference cells.

    gain if mean > max(z * sd_hold / sqrt(n), floor), loss if below the negative; floor is the q
    quantile of |held-out reference mean| over units. Our rule, not part of the tools."""
    sd = hold.std(axis=0)
    floor = float(np.quantile(np.abs(hold.mean(axis=0)), q))
    cut = np.maximum(z * sd / np.sqrt(max(n, 1)), floor)
    return np.where(mean > cut, 3, np.where(mean < -cut, 1, 2))


def clone_profiles(Z, cells, clone, hold):
    """Per (specimen, clone): mean of cell profiles Z (cells x units) and calibrated cn_state."""
    rows = []
    H = Z[hold]
    keep = clone != 'NA'
    spec = cells.specimen.to_numpy()
    for (s, c), idx in pd.Series(np.flatnonzero(keep)).groupby([spec[keep], clone[keep]]):
        m = np.asarray(Z[idx.to_numpy()].mean(axis=0)).ravel()
        st = calibrated_states(m, len(idx), H)
        rows.append(pd.DataFrame(dict(specimen=s, clone=c, bin=np.arange(len(m)), log_ratio=m, cn_state=st)))
    return pd.concat(rows, ignore_index=True)


def differences(profiles, bins, normal_like, min_units=2):
    """Runs of >= min_units consecutive same-chromosome units where two tumour clones differ in cn_state."""
    chrom = bins.set_index('bin').chromosome
    rows = []
    for spec, g in profiles.groupby('specimen'):
        tum = [c for c in g.clone.unique() if c not in normal_like]
        wide_r = g.pivot(index='bin', columns='clone', values='log_ratio')
        wide_s = g.pivot(index='bin', columns='clone', values='cn_state')
        for i, a in enumerate(tum):
            for b in tum[i + 1:]:
                diff = (wide_s[a] != wide_s[b]).to_numpy()
                ub = wide_s.index.to_numpy()
                ch = chrom.reindex(ub).to_numpy()
                start = None
                for k in range(len(ub) + 1):
                    run_on = k < len(ub) and diff[k] and (start is None or ch[k] == ch[start])
                    if run_on and start is None:
                        start = k
                    elif not run_on and start is not None:
                        if k - start >= min_units:
                            sel = ub[start:k]
                            d = float((wide_r.loc[sel, a] - wide_r.loc[sel, b]).mean())
                            bb = bins.set_index('bin').loc[sel]
                            rows.append(dict(specimen=spec, clone_a=a, clone_b=b, chromosome=ch[start],
                                             start=int(bb.start.min()), end=int(bb.end.max()),
                                             n_genes=int(bb.n_genes.sum()), first_bin=int(sel[0]),
                                             last_bin=int(sel[-1]), n_test_bins=len(sel),
                                             direction='higher' if d > 0 else 'lower', mean_diff=d))
                        start = k if (k < len(ub) and diff[k]) else None
    return pd.DataFrame(rows, columns=['specimen', 'clone_a', 'clone_b', 'chromosome', 'start', 'end', 'n_genes',
                                       'first_bin', 'last_bin', 'n_test_bins', 'direction', 'mean_diff'])


def write_layout(out, cells, clone, malignant, bins, genes_df, profiles, normal_like, hold):
    """Write cells.tsv (common), cells.tsv.gz, clones.tsv, bins.tsv, genes.tsv.gz, profiles.tsv.gz,
    arms.tsv, differences.tsv and reference_check.tsv in the Haplodose layout."""
    out = Path(out)
    clone = np.asarray(clone, dtype=object)
    malignant = np.asarray(malignant, dtype=bool)
    write_table(common_cells(cells, clone, malignant), out / 'cells.tsv')
    full = cells.assign(clone=clone, malignant_call=malignant)
    full.to_csv(guard(out / 'cells.tsv.gz'), sep='\t', index=False, compression='gzip')
    bins = bins.copy()
    bins['arm'] = arm_of(bins.chromosome, bins.start)
    write_table(bins, out / 'bins.tsv')
    genes_df.to_csv(guard(out / 'genes.tsv.gz'), sep='\t', index=False, compression='gzip')
    profiles = profiles.assign(log_dose=profiles.log_ratio)[['specimen', 'clone', 'bin', 'log_ratio', 'log_dose', 'cn_state']]
    profiles.to_csv(guard(out / 'profiles.tsv.gz'), sep='\t', index=False, compression='gzip')
    keep = clone != 'NA'
    counts = pd.DataFrame({'specimen': cells.specimen.to_numpy()[keep], 'clone': clone[keep]}).value_counts()
    clones = counts.rename('n_cells_map').reset_index()
    clones['normal_like'] = clones.clone.isin(normal_like)
    write_table(clones[['specimen', 'clone', 'normal_like', 'n_cells_map']], out / 'clones.tsv')
    arm = bins.set_index('bin').arm
    arms = (profiles.assign(arm=profiles.bin.map(arm)).groupby(['specimen', 'clone', 'arm'], sort=False)
            .agg(n_bins=('bin', 'size'), median_log_ratio=('log_ratio', 'median'),
                 mean_log_dose=('log_dose', 'mean')).reset_index())
    write_table(arms, out / 'arms.tsv')
    write_table(differences(profiles, bins, set(normal_like)), out / 'differences.tsv')
    ref = pd.DataFrame({'specimen': cells.specimen.to_numpy()[hold], 'm': malignant[hold]}).groupby('specimen').m
    write_table(pd.DataFrame({'n_tested': ref.size(), 'n_called_aneuploid': ref.sum(),
                              'fraction_called_aneuploid': ref.mean()}).reset_index(), out / 'reference_check.tsv')
