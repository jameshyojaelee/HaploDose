"""Audit corrected original-BAM counts before DNA reconstruction."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record

AUTOSOMES={f'chr{i}' for i in range(1,23)}


def read_counts(path,cells,bins):
    rows={value:i for i,value in enumerate(cells)}
    cols={(str(row.chromosome),int(row.bin)):i for i,row in enumerate(bins.itertuples())}
    matrix=np.zeros((len(cells),len(bins)),np.int64);seen=np.zeros_like(matrix,bool)
    totals=Counter();records=0
    with gzip.open(path,'rt') as stream:
        for rg,chrom,bin_id,number in csv.reader(stream,delimiter='\t'):
            assert rg in rows and bin_id.isdigit() and number.isdigit()
            key=(chrom,int(bin_id));value=int(number);assert value>0
            totals[(rg,chrom)]+=value;records+=1
            if chrom in AUTOSOMES:
                assert key in cols,key
                i,j=rows[rg],cols[key];assert not seen[i,j]
                seen[i,j]=True;matrix[i,j]=value
    return matrix,totals,records


def boundary_flow(old,new):
    delta=new-old
    cumulative=np.cumsum(delta,axis=1)
    assert np.all(cumulative[:,-1]==0),'Chromosome totals changed'
    flow=np.column_stack([np.zeros(len(old),np.int64),cumulative[:,:-1]])
    assert np.all(flow>=0),'Correction is not exclusively leftward'
    assert np.all(flow<=old),'More boundary starts than original source-bin reads'
    incoming=np.column_stack([flow[:,1:],np.zeros(len(old),np.int64)])
    np.testing.assert_array_equal(old-flow+incoming,new)
    return flow


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    root=Path('/gpfs/commons/home/jameslee/Haplodose')
    repair=root/'dev/andor_dna_coordinates_20261108'
    record=repair/'rebuild_23009248';newdir=repair/'corrected_counts_v1'
    assert (record/'COMPLETE').is_file() and 'exit\t0\n' in (record/'job.tsv').read_text()
    assert (newdir/'bins/SNU-601.md5').read_text().strip()=='8c3cf81a44cfd9640034f706fb45671e'
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',rebuild=record/'job.tsv',
                rebuild_hashes=record/'output_sha256.txt',MD5=newdir/'bins/SNU-601.md5',
                old=root/'data/external/andor2020_scdna/bins_original_bam/SNU-601.tsv.gz',
                new=newdir/'bins/SNU-601.tsv.gz',annotation=newdir/'bin_annotation_500kb.tsv',
                old_header=root/'data/external/andor2020_scdna/bins_original_bam/SNU-601.rg.txt',
                new_header=newdir/'bins/SNU-601.rg.txt')
    run=start_run(args.out,inputs,{},dict(scope='Count identity/coordinate flow; no biological validation'),mode='andor_count_repair_audit')
    try:
        for row in inputs['rebuild_hashes'].read_text().splitlines():
            digest,path=row.split(maxsplit=1)
            assert file_record(path)['sha256']==digest
        assert inputs['old_header'].read_bytes()==inputs['new_header'].read_bytes()
        cells=[]
        for line in inputs['new_header'].read_text().splitlines():
            tags=dict(value.split(':',1) for value in line.split('\t')[1:])
            assert tags['ID'] and tags['ID'] not in cells
            cells.append(tags['ID'])
        cells=np.array(sorted(cells))
        ann=pd.read_csv(inputs['annotation'],sep='\t')
        parts=[]
        for chrom,group in ann.groupby('chromosome',sort=False):
            assert group.bin.tolist()==list(range(len(group)))
            # Include one terminal bin to detect historical boundary overflow.
            parts.append(pd.DataFrame(dict(chromosome=chrom,bin=np.arange(len(group)+1))))
        bins=pd.concat(parts,ignore_index=True)
        old,old_totals,nold=read_counts(inputs['old'],cells,bins)
        new,new_totals,nnew=read_counts(inputs['new'],cells,bins)
        assert old_totals==new_totals,'RG/chromosome fragment totals changed'
        np.testing.assert_array_equal(boundary_flow(np.array([[2,3,2]]),np.array([[3,3,1]])),[[0,1,1]])
        flows=np.zeros_like(old);chromosomes=[]
        for chrom,group in bins.groupby('chromosome',sort=False):
            indices=group.index.to_numpy();flow=boundary_flow(old[:,indices],new[:,indices])
            flows[:,indices]=flow
            chromosomes.append(dict(chromosome=chrom,fragments=int(old[:,indices].sum()),
                                     relocated_boundary_fragments=int(flow.sum()),
                                     affected_read_groups=int(np.sum(np.any(flow>0,axis=1)))))
        np.savez_compressed(args.out/'counts.npz',old=old,new=new,flow=flows,cells=cells)
        bins.to_csv(args.out/'bins.tsv',sep='\t',index=False)
        # Resampled and original-depth fits must read identical GC/mask values,
        # without a potentially different CSV float round trip.
        shutil.copy2(inputs['annotation'],args.out/'annotation.tsv')
        (args.out/'read_groups.txt').write_bytes(inputs['new_header'].read_bytes())
        write_json(dict(read_groups=len(cells),autosomal_fragments=int(old.sum()),
                        all_contig_fragments=sum(old_totals.values()),old_rows=nold,new_rows=nnew,
                        relocated_boundary_fragments=int(flows.sum()),
                        affected_read_groups=int(np.sum(np.any(flows>0,axis=1))),
                        per_chromosome=chromosomes,genome_wide_count_repair_verified=True,
                        DNA_clone_impact='not yet reconstructed',biological_specimens=1),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
