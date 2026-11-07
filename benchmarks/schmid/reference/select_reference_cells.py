"""Select the external normal reference cells for claim 5 (docs/STATUS.md item 5, amendment f) and save their IDs.

COLO320: gut cell atlas Full_obj_raw_counts_nosoupx_v2.h5ad; Diagnosis == "Healthy adult", Region == "LargeInt",
  category == "Epithelial", atlas row order; then R `set.seed(1); sample(1:n, 5000)` (computed in R 4.1.0, env
  schmid_r41, by the calling sbatch and passed in as a file of 1-based indices), as in the authors'
  prepare_COLO320_rep1_withRef.R. REC and APD excluded; no rectum sensitivity.
MCF7: Tabula Sapiens TS_Mammary.h5ad; method == "10X", compartment == "epithelial"; all cells, group =
  cell_ontology_class with spaces -> "_", as in prepare_MCF7_dataset_withRef.R.
Only obs annotations are read here (no counts). Usage: select_reference_cells.py OUT R_INDEX_FILE
"""
import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

D = Path('/gpfs/commons/home/jameslee/Haplodose/data/external')
GUT = D / 'gutcellatlas_elmentaite2021/Full_obj_raw_counts_nosoupx_v2.h5ad'
TS = D / 'tabula_sapiens_v1_figshare/TS_Mammary.h5ad'


def decode(v):
    return np.array([x.decode('utf8') if isinstance(x, bytes) else str(x) for x in v])


def col(g, k):
    o = g[k]
    if isinstance(o, h5py.Group) and 'categories' in o:
        cats, codes = decode(o['categories'][:]), o['codes'][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], 'NA')
    if o.attrs.get('categories') is not None:
        cats, codes = decode(g.file[o.attrs['categories']][:]), o[:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], 'NA')
    v = o[:]
    return decode(v) if v.dtype.kind in 'SO' else v


def obs(path, keys):
    with h5py.File(path) as h:
        o = h['obs']
        idx = o.attrs.get('_index', '_index')
        return pd.DataFrame({k: col(o, k) for k in keys}, index=pd.Index(col(o, idx), name='cell_id'))


def md5(p):
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


def main(out, r_index_file):
    out = Path(out)
    g = obs(GUT, ['Diagnosis', 'Region', 'category', 'Age_group', '10X', 'Sample name'])
    g['atlas_row'] = np.arange(len(g)) + 1
    sub = g[g.Diagnosis.eq('Healthy adult') & g.Region.eq('LargeInt') & g.category.eq('Epithelial')]
    idx = np.loadtxt(r_index_file, dtype=int)
    assert len(idx) == 5000 and len(set(idx)) == 5000 and idx.min() >= 1 and idx.max() <= len(sub)
    sel = sub.iloc[idx - 1].assign(r_index=idx, group='Epithelial')
    sel.reset_index().to_csv(out / 'COLO320_reference_cells.tsv', sep='\t', index=False)

    t = obs(TS, ['method', 'compartment', 'cell_ontology_class', 'donor'])
    t['atlas_row'] = np.arange(len(t)) + 1
    ts = t[t.method.eq('10X') & t.compartment.eq('epithelial')].copy()
    ts['group'] = ts.cell_ontology_class.str.replace(' ', '_', regex=False)
    ts.reset_index().to_csv(out / 'MCF7_reference_cells.tsv', sep='\t', index=False)

    summary = dict(
        colo320=dict(source=str(GUT), md5=md5(GUT), n_eligible=len(sub), n_selected=len(sel),
                     rule='Diagnosis=="Healthy adult" & Region=="LargeInt" & category=="Epithelial", atlas order; '
                          'R 4.1.0 set.seed(1); sample(1:n, 5000)',
                     chemistry=sel['10X'].value_counts().to_dict(), samples=int(sel['Sample name'].nunique()),
                     age_group=sel.Age_group.value_counts().to_dict()),
        mcf7=dict(source=str(TS), md5=md5(TS), n_selected=len(ts), rule='method=="10X" & compartment=="epithelial"',
                  groups=ts.group.value_counts().to_dict(), donors=ts.donor.unique().tolist()))
    (out / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main(*sys.argv[1:3])
