"""Read-only: gene identifiers in the normal-section h5ads vs the Space Ranger panel."""
import h5py
import numpy as np
import pandas as pd
import anndata as ad
with h5py.File('/gpfs/commons/groups/innovation/share/Spaceranger/novaseq/U0396/outs/segmented_outputs/filtered_feature_cell_matrix.h5', 'r') as h:
    ids = pd.Index(h['matrix/features/id'].asstr()[:]); names = pd.Index(h['matrix/features/name'].asstr()[:])
print('panel ids', ids[:3].tolist(), 'unique', ids.is_unique, 'names unique', names.is_unique, 'dup names', names[names.duplicated()].tolist()[:20])
for s in ['JGCRC4NT', 'JGCRC14NT', 'P3_NAT', 'P5_NAT']:
    a = ad.read_h5ad(f'/gpfs/commons/groups/innovation/share/h5ad/processed_h5ad_090126/{s}_processed.h5ad', backed='r')
    v = a.var
    print(s, v.head(3).to_dict(), 'var_names in panel names', np.isin(a.var_names, names).mean(),
          'panel names in var_names', names.isin(a.var_names).mean(), 'gene_ids in panel ids', np.isin(v.gene_ids.astype(str), ids).mean(),
          'same order as panel', len(v) == len(names) and bool(np.array_equal(a.var_names, names)))
    extra = a.var_names[~np.isin(a.var_names, names)]
    print('  not in panel (first 15):', extra[:15].tolist())
    a.file.close()
