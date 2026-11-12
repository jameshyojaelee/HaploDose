"""DNA-only reconstruction and label-invariant comparisons of uncertain groups."""
import gzip
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score

from haplodose.provenance import file_record


def dna_module():
    path=Path(__file__).resolve().parents[1]/'andor/scdna_cn.py'
    spec=importlib.util.spec_from_file_location('frozen_dna_reconstruction',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def verify_run(directory):
    run=json.loads((directory/'run_manifest.json').read_text())
    assert run['status']=='complete'
    for group in ('input','output'):
        for record in run[group]['files'].values():
            assert file_record(record['path'])['sha256']==record['sha256']
    inventory=directory/'output_files.json'
    if inventory.exists():
        for record in json.loads(inventory.read_text()).values():
            assert file_record(record['path'])['sha256']==record['sha256']
    return run


def write_inventory(directory):
    from haplodose.io import write_json
    records={str(path.relative_to(directory)):file_record(path) for path in sorted(directory.rglob('*'))
             if path.is_file() and path not in (directory/'run_manifest.json',directory/'output_files.json')}
    write_json(records,directory/'output_files.json')


def fit(module,bins,destination):
    module.main(['--bins-dir',str(bins),'--lines','SNU-601','--out',str(destination),
                 '--min-fragments','50000','--max-mapd','0.45','--t-split','5',
                 '--min-diff-mb','10','--min-clone','10'])
    return destination/'SNU-601'


def assigned(frame):
    return ~frame['clone'].isin(['excluded','unassigned'])


def compare(anchor,candidate):
    a=pd.read_csv(anchor/'cells.tsv',sep='\t',dtype={'rg':str}).set_index('rg')
    b=pd.read_csv(candidate/'cells.tsv',sep='\t',dtype={'rg':str}).set_index('rg')
    assert a.index.is_unique and b.index.is_unique
    common=a.index.intersection(b.index)
    both=common[assigned(a.loc[common]).to_numpy()&assigned(b.loc[common]).to_numpy()]
    matched=0
    if len(both):
        table=pd.crosstab(a.loc[both,'clone'],b.loc[both,'clone']).to_numpy()
        row,col=linear_sum_assignment(-table);matched=int(table[row,col].sum())
    aq=json.loads((anchor/'qc.json').read_text());bq=json.loads((candidate/'qc.json').read_text())
    # These are outputs generated locally by the recorded DNA writer, which
    # stores string cell IDs as object arrays. Do not apply this loader to
    # arbitrary external NPZs. Verify ID membership against the QC cell table;
    # each channel can retain its own row order and is aligned below by ID.
    with np.load(anchor/'cn.npz',allow_pickle=True) as archive:
        az={name:archive[name] for name in ('cells','bins','cn','seg_ratio')}
    with np.load(candidate/'cn.npz',allow_pickle=True) as archive:
        bz={name:archive[name] for name in ('cells','bins','cn','seg_ratio')}
    assert pd.Index(az['cells']).is_unique and pd.Index(bz['cells']).is_unique
    np.testing.assert_array_equal(np.sort(az['cells']),np.sort(a.index[a.pass_mapd].to_numpy()))
    np.testing.assert_array_equal(np.sort(bz['cells']),np.sort(b.index[b.pass_mapd].to_numpy()))
    np.testing.assert_array_equal(az['bins'],bz['bins'])
    ids=pd.Index(az['cells']).intersection(pd.Index(bz['cells']))
    ia=pd.Index(az['cells']).get_indexer(ids);ib=pd.Index(bz['cells']).get_indexer(ids)
    relative=np.abs(az['seg_ratio'][ia].astype(float)-bz['seg_ratio'][ib].astype(float))
    integer=az['cn'][ia]!=bz['cn'][ib]
    arms=lambda d:sorted(set(pd.read_csv(d/'clone_differences.tsv',sep='\t').arm))
    aa,ba=arms(anchor),arms(candidate)
    return dict(anchor_cells=len(a),candidate_cells=len(b),common_cells=len(common),
                anchor_assigned=int(assigned(a).sum()),candidate_assigned=int(assigned(b).sum()),
                both_assigned=len(both),one_to_one_matched_assigned=matched,
                matched_assigned_fraction_of_anchor_cells=matched/len(a),
                ARI_on_both_assigned=float(adjusted_rand_score(a.loc[both,'clone'],b.loc[both,'clone'])) if len(both) else None,
                anchor_pass_depth=int(a.pass_depth.sum()),candidate_pass_depth=int(b.pass_depth.sum()),
                anchor_pass_MAPD=int(a.pass_mapd.sum()),candidate_pass_MAPD=int(b.pass_mapd.sum()),
                anchor_clones=len(aq['clones']),candidate_clones=len(bq['clones']),
                anchor_ploidy_scale=aq['ploidy_scale'],candidate_ploidy_scale=bq['ploidy_scale'],
                DNA_profile_cells_compared=len(ids),
                mean_absolute_relative_profile_change=float(relative.mean()) if len(ids) else None,
                maximum_absolute_relative_profile_change=float(relative.max()) if len(ids) else None,
                fraction_integer_CN_entries_changed=float(integer.mean()) if len(ids) else None,
                anchor_inferred_difference_arms=aa,candidate_inferred_difference_arms=ba,
                shared_inferred_difference_arms=sorted(set(aa)&set(ba)),
                interpretation='Inferred DNA measurement agreement, not genomic truth or RNA membership accuracy')


def write_counts(path,matrix,cells,bins):
    chrom=bins.chromosome.to_numpy();bin_id=bins.bin.to_numpy()
    with gzip.open(path,'wt',compresslevel=1) as stream:
        for i,rg in enumerate(cells):
            keep=np.flatnonzero(matrix[i])
            stream.writelines(f'{rg}\t{chrom[j]}\t{bin_id[j]}\t{int(matrix[i,j])}\n' for j in keep)
