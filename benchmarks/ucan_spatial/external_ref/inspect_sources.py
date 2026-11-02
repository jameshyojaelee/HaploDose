"""Read-only inspection of the normal-section sources: h5ad obs, matrix barcodes, gene panels."""
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

H5AD = Path('/gpfs/commons/groups/innovation/share/h5ad/processed_h5ad_090126')
B2C = Path('/gpfs/commons/groups/innovation/share/10X_VisiumHD_CRC/adata/annotated_b2c_adatas')
SEG = {
    'JGCRC4NT': Path('/gpfs/commons/groups/innovation/share/Spaceranger/novaseq/JGCRC4NT/outs/segmented_outputs'),
    'JGCRC14NT': Path('/gpfs/commons/groups/innovation/share/Spaceranger/ultima_v4/JGCRC14NT/outs/segmented_outputs'),
}
pd.set_option('display.width', 200)


def show_adata(path):
    a = ad.read_h5ad(path, backed='r')
    print('##', path, a.shape)
    print('obs columns:', a.obs.columns.tolist())
    print('obs_names head:', a.obs_names[:5].tolist())
    print('layers:', list(a.layers.keys()) if hasattr(a, 'layers') else None, 'raw:', a.raw is not None)
    print('obsm:', list(a.obsm.keys()))
    for c in ['cell_type', 'rctd_top_first_type', 'annotation_main']:
        if c in a.obs:
            print(c, a.obs[c].value_counts(dropna=False).head(40).to_dict())
    if 'cell_type' in a.obs and 'annotation_main' in a.obs:
        print(pd.crosstab(a.obs.cell_type, a.obs.annotation_main).to_string())
    for c in a.obs.columns:
        if any(k in c.lower() for k in ['px', 'py', 'x', 'y', 'count', 'umi', 'object', 'bin', 'id']):
            print('  col', c, a.obs[c].dtype, a.obs[c].iloc[:3].tolist())
    X = a.X[:5]
    X = X.toarray() if hasattr(X, 'toarray') else np.asarray(X)
    print('X dtype', X.dtype, 'first rows sum', X.sum(axis=1), 'integer', np.allclose(X, np.round(X)))
    print('var head', a.var_names[:5].tolist(), a.var.columns.tolist())
    a.file.close()


for s in ['JGCRC4NT', 'JGCRC14NT', 'P3_NAT', 'P5_NAT']:
    show_adata(H5AD / f'{s}_processed.h5ad')
for p in sorted(B2C.glob('*NAT*')) + [B2C / 'JGCRC14NT_annotated.h5ad']:
    show_adata(p)

for s, d in SEG.items():
    with h5py.File(d / 'filtered_feature_cell_matrix.h5', 'r') as h:
        bc = h['matrix/barcodes'].asstr()[:]
        genes = h['matrix/features/name'].asstr()[:]
        print('##', s, 'segmented cells', len(bc), bc[:3], 'genes', len(genes), 'IGF2' in genes)
        print(' attrs', {k: str(v) for k, v in h.attrs.items()})
    print(' scale', json.loads((d / 'spatial/scalefactors_json.json').read_text()))
with h5py.File('/gpfs/commons/home/jameslee/CRC/results/cnv_caller_benchmark/U0396/input/counts.h5', 'r') as h:
    g = h['genes'].asstr()[:]
    print('## U0396 input genes', len(g), g[:5], 'cells', len(h['cell_ids']))
with h5py.File('/gpfs/commons/groups/innovation/share/Spaceranger/novaseq/U0396/outs/segmented_outputs/filtered_feature_cell_matrix.h5', 'r') as h:
    print(' U0396 seg attrs', {k: str(v) for k, v in h.attrs.items()})
for s in ['P3', 'P5']:
    d = Path(f'/gpfs/commons/groups/innovation/share/10X_VisiumHD_CRC/data/{s}_NAT/binned_outputs/square_002um')
    with h5py.File(d / 'filtered_feature_bc_matrix.h5', 'r') as h:
        print('##', s, '2um bins', len(h['matrix/barcodes']), h['matrix/barcodes'].asstr()[:3],
              'genes', len(h['matrix/features/name']), {k: str(v) for k, v in h.attrs.items()})
    print(' scale', json.loads((d / 'spatial/scalefactors_json.json').read_text()))
m = pd.read_csv('/gpfs/commons/home/jameslee/CRC/results/cnv_caller_benchmark/U0396/input/cells.tsv', sep='\t')
print(m.role.value_counts().to_dict(), m.groupby(['ref_group', 'ref_split']).size().to_dict())
print(m[m.role.eq('reference')].total_umi_included.describe())
print(m[m.role.eq('tumour_epithelium')].total_umi_included.describe())
