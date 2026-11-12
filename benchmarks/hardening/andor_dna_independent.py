"""Independent integer/rational references for locally generated DNA artifacts."""
from collections import Counter
import csv
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path

import networkx as nx
import numpy as np


def choose2(n):
    return n*(n-1)//2


def partition_reference(left,right):
    assert len(left)==len(right)
    if not left:
        return dict(matched=0,ARI=None)
    table=Counter(zip(left,right));rows=Counter(left);cols=Counter(right)
    graph=nx.Graph()
    for (a,b),count in table.items():
        graph.add_edge(('left',a),('right',b),weight=int(count))
    pairs=nx.max_weight_matching(graph,weight='weight')
    matched=sum(graph[a][b]['weight'] for a,b in pairs)
    n=len(left);total=choose2(n)
    observed=sum(choose2(value) for value in table.values())
    row_pairs=sum(choose2(value) for value in rows.values())
    col_pairs=sum(choose2(value) for value in cols.values())
    if total:
        expectation=Fraction(row_pairs*col_pairs,total)
        maximum=Fraction(row_pairs+col_pairs,2)
        ari=(Fraction(observed)-expectation)/(maximum-expectation) if maximum!=expectation else Fraction(1)
    else:
        ari=Fraction(1)
    return dict(matched=int(matched),ARI=float(ari),ARI_fraction=str(ari))


def read_cells(directory):
    result={}
    with (directory/'cells.tsv').open() as stream:
        for row in csv.DictReader(stream,delimiter='\t'):
            assert row['rg'] and row['rg'] not in result
            for field in ('pass_depth','pass_mapd'):
                assert row[field] in ('True','False')
                row[field]=row[field]=='True'
            result[row['rg']]=row
    return result


def arm_reference(directory):
    qc=json.loads((directory/'qc.json').read_text());names=list(qc['clones'])
    with (directory/'clone_bins.tsv').open() as stream:
        bins=list(csv.DictReader(stream,delimiter='\t'))
    expected={}
    for left,right in itertools.combinations(names,2):
        for arm in dict.fromkeys(row['arm'] for row in bins):
            rows=[row for row in bins if row['arm']==arm]
            supported=0;different=0;weighted=Fraction(0)
            for row in rows:
                a,b=row[left],row[right]
                if not a or not b or not math.isfinite(float(a)) or not math.isfinite(float(b)):
                    continue
                length=int(row['end'])-int(row['start']);assert length>0
                delta=Fraction(a)-Fraction(b);supported+=length
                if abs(delta)>=1:
                    different+=length;weighted+=delta*length
            if supported and 10*different>=3*supported:
                expected[(left,right,arm)]=(float(Fraction(different,supported)),float(weighted/different))
    with (directory/'clone_differences.tsv').open() as stream:
        actual=list(csv.DictReader(stream,delimiter='\t'))
    assert len(actual)==len(expected)
    seen=set()
    for row in actual:
        key=(row['clone_a'],row['clone_b'],row['arm']);assert key in expected and key not in seen
        seen.add(key)
        assert float(row['frac_bins_diff'])==expected[key][0]
        assert float(row['mean_diff'])==expected[key][1]
    return sorted(set(key[2] for key in expected))


def audit_comparison(anchor,candidate,reported):
    a,b=read_cells(anchor),read_cells(candidate)
    is_assigned=lambda row:row['clone'] not in ('excluded','unassigned')
    ids=[key for key in a if key in b]
    both=[key for key in ids if is_assigned(a[key]) and is_assigned(b[key])]
    reference=partition_reference([a[key]['clone'] for key in both],[b[key]['clone'] for key in both])
    aq=json.loads((anchor/'qc.json').read_text());bq=json.loads((candidate/'qc.json').read_text())
    for frame,qc in ((a,aq),(b,bq)):
        assert dict(Counter(row['clone'] for row in frame.values() if is_assigned(row)))==qc['clones']
        assert sum(row['clone']=='unassigned' for row in frame.values())==qc['n_unassigned']
        assert sum(row['pass_depth'] for row in frame.values())==qc['n_pass_depth']
        assert sum(row['pass_mapd'] for row in frame.values())==qc['n_pass_mapd']
    aa,ba=arm_reference(anchor),arm_reference(candidate)
    expected=dict(anchor_cells=len(a),candidate_cells=len(b),common_cells=len(ids),
                  anchor_assigned=sum(is_assigned(row) for row in a.values()),
                  candidate_assigned=sum(is_assigned(row) for row in b.values()),
                  both_assigned=len(both),one_to_one_matched_assigned=reference['matched'],
                  matched_assigned_fraction_of_anchor_cells=reference['matched']/len(a),
                  ARI_on_both_assigned=reference['ARI'],
                  anchor_pass_depth=sum(row['pass_depth'] for row in a.values()),
                  candidate_pass_depth=sum(row['pass_depth'] for row in b.values()),
                  anchor_pass_MAPD=sum(row['pass_mapd'] for row in a.values()),
                  candidate_pass_MAPD=sum(row['pass_mapd'] for row in b.values()),
                  anchor_clones=len(aq['clones']),candidate_clones=len(bq['clones']),
                  anchor_ploidy_scale=aq['ploidy_scale'],candidate_ploidy_scale=bq['ploidy_scale'],
                  anchor_inferred_difference_arms=aa,candidate_inferred_difference_arms=ba,
                  shared_inferred_difference_arms=sorted(set(aa)&set(ba)))
    for key,value in expected.items():
        assert reported[key]==value,(key,reported[key],value)
    # Files are generated locally by the hash-verified DNA writer; legacy IDs
    # are object arrays. Never use this loader for arbitrary external NPZs.
    with np.load(anchor/'cn.npz',allow_pickle=True) as archive:
        az={key:archive[key] for key in ('cells','bins','cn','seg_ratio')}
    with np.load(candidate/'cn.npz',allow_pickle=True) as archive:
        bz={key:archive[key] for key in ('cells','bins','cn','seg_ratio')}
    assert len(set(az['cells']))==len(az['cells']) and len(set(bz['cells']))==len(bz['cells'])
    assert set(az['cells'])=={key for key,row in a.items() if row['pass_mapd']}
    assert set(bz['cells'])=={key for key,row in b.items() if row['pass_mapd']}
    np.testing.assert_array_equal(az['bins'],bz['bins'])
    ai={value:i for i,value in enumerate(az['cells'])};bi={value:i for i,value in enumerate(bz['cells'])}
    shared=sorted(set(ai)&set(bi));assert len(shared)==reported['DNA_profile_cells_compared']
    sums=[];source_scales=[];maximum=0.;changed=0;entries=0
    for cell in shared:
        x=az['seg_ratio'][ai[cell]].astype(float);y=bz['seg_ratio'][bi[cell]].astype(float)
        assert np.isfinite(x).all() and np.isfinite(y).all()
        differences=np.abs(x-y)
        sums.append(math.fsum(differences));source_scales.append(math.fsum(np.abs(x)+np.abs(y)))
        maximum=max(maximum,float(max(differences)))
        changed+=sum(int(left!=right) for left,right in zip(az['cn'][ai[cell]],bz['cn'][bi[cell]]))
        entries+=len(x)
    bound=0.
    if entries:
        exact_sum=math.fsum(sums);mean=exact_sum/entries
        epsilon=np.finfo(float).eps/2
        gamma=(entries-1)*epsilon/(1-(entries-1)*epsilon)
        # Nonnegative reduction bound, subtraction rounding, and final division.
        bound=gamma*mean+epsilon*math.fsum(source_scales)/entries+2*math.ulp(mean)
        assert abs(reported['mean_absolute_relative_profile_change']-mean)<=bound
        assert reported['maximum_absolute_relative_profile_change']==maximum
        assert reported['fraction_integer_CN_entries_changed']==changed/entries
    else:
        assert all(reported[key] is None for key in ('mean_absolute_relative_profile_change',
                    'maximum_absolute_relative_profile_change','fraction_integer_CN_entries_changed'))
    return dict(cell_counts_and_matching_exact=True,ARI_fraction=reference.get('ARI_fraction'),
                arm_contrasts_reconstructed=True,profile_entries=entries,
                profile_mean_rounding_bound=bound)
