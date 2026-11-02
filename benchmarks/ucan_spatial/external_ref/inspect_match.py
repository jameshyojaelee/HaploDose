"""Read-only: JGCRC4NT region labels, and h5ad obs vs Space Ranger segmented barcode matching."""
import anndata as ad
import h5py
import numpy as np
import pandas as pd

H5AD = '/gpfs/commons/groups/innovation/share/h5ad/processed_h5ad_090126/{}_processed.h5ad'
SEG = {'JGCRC4NT': '/gpfs/commons/groups/innovation/share/Spaceranger/novaseq/JGCRC4NT/outs/segmented_outputs',
       'JGCRC14NT': '/gpfs/commons/groups/innovation/share/Spaceranger/ultima_v4/JGCRC14NT/outs/segmented_outputs',
       'U0396': '/gpfs/commons/groups/innovation/share/Spaceranger/novaseq/U0396/outs/segmented_outputs'}
pd.set_option('display.width', 250)
a = ad.read_h5ad(H5AD.format('JGCRC4NT'), backed='r')
for c in ['manual_clustering', 'labels_joint_source']:
    print(c, a.obs[c].value_counts().head(30).to_dict())
print(pd.crosstab(a.obs.cell_type, a.obs.manual_clustering).to_string())
a.file.close()
for s, d in SEG.items():
    a = ad.read_h5ad(H5AD.format(s), backed='r')
    obs = a.obs
    with h5py.File(d + '/filtered_feature_cell_matrix.h5', 'r') as h:
        bc = pd.Index(h['matrix/barcodes'].asstr()[:])
        ip = h['matrix/indptr'][:]
        data = h['matrix/data']
        tot = np.add.reduceat(data[:], ip[:-1]) if len(data) else np.zeros(len(bc))
        tot[np.diff(ip) == 0] = 0
    num = pd.Index(bc.str.replace('cellid_', '').str.replace('-1', '').astype(int))
    oid = obs.object_id.astype(int).to_numpy()
    hit = num.get_indexer(oid)
    m = hit >= 0
    print(f'## {s}: h5ad cells {len(obs)}, segmented {len(bc)}, h5ad object_id found in segmented {m.sum()} '
          f'({m.mean():.3f}); segmented found in h5ad {np.isin(num, oid).mean():.3f}')
    if m.any():
        n = obs.nCount_RNA.to_numpy()[m]
        print('   nCount_RNA equal segmented total (all genes):', np.mean(n == tot[hit[m]]),
              'median ratio', np.median(n / np.maximum(tot[hit[m]], 1)))
    a.file.close()
