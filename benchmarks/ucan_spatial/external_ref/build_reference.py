"""External normal colon epithelium reference for the U-CAN Visium HD runs (40 um tiles).

Sources (read only): /gpfs/commons/groups/innovation/share/h5ad/processed_h5ad_090126/<S>_processed.h5ad,
X = raw integer counts per segmented cell (bin2cell objects; these are NOT the Space Ranger segmented
cells: object ids do not match the segmented matrix barcodes, see barcode_match.tsv).

Cell selection per section:
  - cell_type in NORMAL_EPI (normal epithelial labels; all Tumor_* labels excluded)
  - rctd_top_first_type == 'Epi'
  - annotation_main == 'Mucosa' where the column exists (JGCRC14NT, P3_NAT, P5_NAT; removes the
    JGCRC14NT Tumor regions). JGCRC4NT has no annotation_main: keep only epithelial manual_clustering
    clusters (EpiT_Epi_*) in which normal-epithelial labels outnumber Tumor_* labels (stand-in for the
    missing Mucosa region).
Genes: h5ad genes outside the 18,085-gene filtered panel (non-included probes in JGCRC4NT/14NT) dropped.
Tiling as CRC src/43_cnv_caller/spatial.py prepare: counts summed over cells in the same 40 um square
(floor(position_um / 40)) and the same compartment (here one compartment, normal_epithelium); tile
x, y = mean cell position. Position in um = array_col / array_row (2 um bin grid) x 2.
Tiles need >= 300 UMIs (genes aligned to the U-CAN 18,085-gene panel). Per section at most
REF_NORMAL_EPI = 1500 tiles (CRC CNV_infercnv cap for normal epithelium), split fit/holdout 70/30.
Usage: build_reference.py OUTDIR UCAN_INPUT_COUNTS_H5 (gene order template)
"""
import json
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

SEED = 20261101
TILE_UM = 40.0
MIN_UMI = 300
CAP = 1500
FIT_FRACTION = 0.7
OFFSET = 1e8   # um; external tiles sit far from any section (Haplodose's spatial graph uses query rows only)
SECTIONS = ['JGCRC4NT', 'JGCRC14NT', 'P3_NAT', 'P5_NAT']
H5AD = Path('/gpfs/commons/groups/innovation/share/h5ad/processed_h5ad_090126')
SEG = {'JGCRC4NT': '/gpfs/commons/groups/innovation/share/Spaceranger/novaseq/JGCRC4NT/outs/segmented_outputs',
       'JGCRC14NT': '/gpfs/commons/groups/innovation/share/Spaceranger/ultima_v4/JGCRC14NT/outs/segmented_outputs'}
NORMAL_EPI = ['Enterocyte', 'Goblet', 'Immature_Goblet', 'Best4', 'Enteroendocrine', 'Stem_TAlike', 'Tuft']

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=False)
with h5py.File(sys.argv[2], 'r') as h:
    panel = pd.Index(h['genes'].asstr()[:])
# Three gene names repeat in the panel (TBCE, HSPA14, TMSB15B). The h5ads use anndata's
# var_names_make_unique (second copy 'TBCE-1'), so build the same unique names in panel order.
seen = pd.Series(panel).groupby(panel).cumcount().to_numpy()
panel_unique = pd.Index([g if k == 0 else f'{g}-{k}' for g, k in zip(panel, seen)])
if not panel_unique.is_unique:
    raise ValueError('panel names not unique after make_unique')
rng = np.random.default_rng(SEED)
blocks, metas, report, match = [], [], [], []
for s in SECTIONS:
    a = ad.read_h5ad(H5AD / f'{s}_processed.h5ad')
    obs = a.obs
    label = obs.cell_type.astype(str)
    tumour_label = label.str.startswith('Tumor_')
    keep = label.isin(NORMAL_EPI) & obs.rctd_top_first_type.astype(str).eq('Epi')
    if 'annotation_main' in obs:
        region_rule = 'annotation_main == Mucosa'
        keep &= obs.annotation_main.astype(str).eq('Mucosa')
    else:
        cl = obs.manual_clustering.astype(str)
        n_norm = label.isin(NORMAL_EPI).groupby(cl).sum()
        n_tum = tumour_label.groupby(cl).sum()
        ok = n_norm.index[(n_norm > n_tum) & n_norm.index.str.startswith('EpiT_Epi')]
        region_rule = 'manual_clustering in ' + ','.join(ok)
        keep &= cl.isin(ok)
    # Barcode check against the Space Ranger segmented matrix (cellid_<9-digit id>-1).
    if s in SEG:
        with h5py.File(SEG[s] + '/filtered_feature_cell_matrix.h5', 'r') as h:
            bc = pd.Index(h['matrix/barcodes'].asstr()[:])
            ip = h['matrix/indptr'][:]
            tot = np.zeros(len(bc))
            nz = np.diff(ip) > 0
            tot[nz] = np.add.reduceat(h['matrix/data'][:], ip[:-1][nz])
        want = 'cellid_' + obs.object_id.astype(int).map('{:09d}'.format) + '-1'
        hit = bc.get_indexer(want)
        found = hit >= 0
        same = found & (obs.nCount_RNA.to_numpy() == np.where(found, tot[hit], -1))
        match.append(dict(section=s, h5ad_cells=len(obs), segmented_cells=len(bc), id_found=int(found.sum()),
                          id_found_fraction=float(found.mean()), id_found_and_same_umi=int(same.sum()),
                          same_umi_fraction=float(same.mean())))
    else:
        match.append(dict(section=s, h5ad_cells=len(obs), segmented_cells=0, id_found=0, id_found_fraction=np.nan,
                          id_found_and_same_umi=0, same_umi_fraction=np.nan))
    sel = np.flatnonzero(keep.to_numpy())
    X = a.X[sel]
    X = sparse.csr_matrix(X) if not sparse.issparse(X) else X.tocsr()
    col = panel_unique.get_indexer(a.var_names)
    n_missing = int((~panel_unique.isin(a.var_names)).sum())
    X = X[:, col >= 0].tocoo()
    X = sparse.csr_matrix((X.data, (X.row, col[col >= 0][X.col])), shape=(len(sel), len(panel)))
    um = obs.iloc[sel][['array_col', 'array_row']].to_numpy(np.float64) * 2.0
    tile = np.floor(um / TILE_UM).astype(np.int64)
    keys = pd.Index([f'{s}_{x}_{y}_normal_epithelium' for x, y in tile])
    names = pd.Index(keys.unique())
    g = names.get_indexer(keys)
    agg = sparse.csr_matrix((np.ones(len(sel)), (g, np.arange(len(sel)))), shape=(len(names), len(sel)))
    T = (agg @ X).tocsr()
    if int(T.sum()) != int(X.sum()):
        raise ValueError('tile aggregation lost counts')
    umi = np.asarray(T.sum(axis=1)).ravel()
    xy = pd.DataFrame(um, columns=['x', 'y']).groupby(g).mean().reindex(range(len(names)))
    n_cells = np.bincount(g, minlength=len(names))
    ok = np.flatnonzero(umi >= MIN_UMI)
    pick = np.sort(rng.choice(ok, size=min(len(ok), CAP), replace=False))
    order = rng.permutation(len(pick))
    split = np.where(order < int(round(FIT_FRACTION * len(pick))), 'fit', 'holdout')
    meta = pd.DataFrame(dict(cell_id=names[pick], compartment_strict='normal_epithelium', IGF2_UMI=0,
                             total_umi_included=umi[pick].astype(np.int64),
                             x=xy.x.to_numpy()[pick] + OFFSET * (1 + SECTIONS.index(s)), y=xy.y.to_numpy()[pick] + OFFSET, n_segmented_cells=n_cells[pick],
                             dataset='VisiumHD', patient=s, sample=s, specimen=s, role='reference',
                             ref_group=f'normal_epi_{s}', ref_split=split))
    blocks.append(T[pick])
    metas.append(meta)
    report.append(dict(section=s, region_rule=region_rule, h5ad_cells=len(obs), selected_cells=len(sel),
                       selected_by_type=label.iloc[sel].value_counts().to_dict(), panel_genes_missing=n_missing,
                       h5ad_genes_not_in_panel=int((col < 0).sum()), tiles=len(names), tiles_ge_300=len(ok),
                       tiles_kept=len(pick), fit=int((split == 'fit').sum()), holdout=int((split == 'holdout').sum()),
                       median_umi_kept=float(np.median(umi[pick])), median_cells_per_tile=float(np.median(n_cells[pick]))))
    print(json.dumps(report[-1]), flush=True)
    del a, X, T
C = sparse.vstack(blocks).tocsr()
C.sort_indices()
meta = pd.concat(metas, ignore_index=True)
with h5py.File(out / 'counts.h5', 'w') as h:
    h.create_dataset('genes', data=panel.to_numpy().astype(h5py.string_dtype()))
    h.create_dataset('cell_ids', data=meta.cell_id.to_numpy(dtype=h5py.string_dtype()))
    grp = h.create_group('X')
    for key, values in [('data', C.data.astype(np.int64)), ('indices', C.indices), ('indptr', C.indptr)]:
        grp.create_dataset(key, data=values, compression='gzip')
meta.to_csv(out / 'cells.tsv', sep='\t', index=False)
pd.DataFrame(match).to_csv(out / 'barcode_match.tsv', sep='\t', index=False)
(out / 'build.json').write_text(json.dumps(dict(seed=SEED, tile_um=TILE_UM, min_umi=MIN_UMI, cap=CAP,
                                                fit_fraction=FIT_FRACTION, normal_epi_labels=NORMAL_EPI,
                                                gene_template=sys.argv[2], sections=report), indent=2) + '\n')
print(pd.DataFrame(match).to_string(index=False))
