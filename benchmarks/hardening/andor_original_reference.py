"""Check the original DNA BAM reference contract without reprocessing its reads."""
import argparse
import json
from pathlib import Path
import subprocess

import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def sequence_dictionary(text):
    result={}
    for line in text.splitlines():
        if line.startswith('@SQ\t'):
            tags=dict(field.split(':',1) for field in line.split('\t')[1:])
            assert tags['SN'] not in result
            result[tags['SN']]=tags
    assert result
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    root=Path('/gpfs/commons/home/jameslee/Haplodose')
    fasta=Path('/gpfs/commons/groups/vickovic_lab/refs/human/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa')
    annotation=root/'dev/andor_dna_coordinates_20261108/corrected_counts_v1/bin_annotation_500kb.tsv'
    sources=root/'dev/andor_dna_coordinates_20261108/rebuild_23009248/bam_sources.tsv'
    run=start_run(args.out,dict(script=Path(__file__),FAI=Path(str(fasta)+'.fai'),
                  annotation=annotation,sources=sources),{},dict(scope='Original DNA reference identity; no RNA or clone accuracy claim'),
                  mode='andor_original_reference')
    try:
        versions=subprocess.run(['samtools','--version'],capture_output=True,text=True,check=True)
        (args.out/'samtools_version.txt').write_text(versions.stdout)
        table=pd.read_csv(sources,sep='\t');row=table.loc[table.line=='SNU-601'].iloc[0]
        command=['samtools','view','--no-PG','-H',row.url]
        header=subprocess.run(command,capture_output=True,text=True,timeout=300)
        (args.out/'header.sam').write_text(header.stdout)
        (args.out/'header.stderr').write_text(header.stderr)
        header.check_returncode()
        observed=sequence_dictionary(header.stdout)
        fai=pd.read_csv(str(fasta)+'.fai',sep='\t',header=None)
        lengths=dict(zip(fai[0],fai[1]));ann=pd.read_csv(annotation,sep='\t')
        records=[]
        for chrom in (f'chr{i}' for i in range(1,23)):
            assert chrom in observed and chrom in lengths
            length=int(lengths[chrom]);assert int(observed[chrom]['LN'])==length
            bins=ann.loc[ann.chromosome==chrom].sort_values('bin')
            assert bins.bin.tolist()==list(range((length+499999)//500000))
            assert bins.start.tolist()==list(range(0,length,500000))
            assert bins.end.tolist()==[min(start+500000,length) for start in range(0,length,500000)]
            records.append(dict(chromosome=chrom,bases=length,BAM_M5=observed[chrom].get('M5')))
        have_md5=all(record['BAM_M5'] for record in records)
        if have_md5:
            command_dict=['samtools','dict',str(fasta)]
            dictionary=subprocess.run(command_dict,capture_output=True,text=True,check=True,timeout=1200)
            (args.out/'local_reference.dict').write_text(dictionary.stdout)
            (args.out/'local_reference.stderr').write_text(dictionary.stderr)
            local=sequence_dictionary(dictionary.stdout)
            for record in records:
                assert record['BAM_M5']==local[record['chromosome']]['M5']
            write_json(file_record(fasta),args.out/'reference_file.json')
        write_json(dict(commands=[command]+([command_dict] if have_md5 else []),
                        autosomal_lengths_and_annotation_intervals_exact=True,
                        autosomal_sequence_MD5_identical=have_md5,chromosomes=records,
                        reference_sequence_identity='verified by BAM/local sequence MD5' if have_md5 else
                        'unverified: BAM header lacks complete autosomal sequence checksums',
                        biological_specimens=1,independent_biological_validation=False),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
