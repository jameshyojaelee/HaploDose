"""Research metadata-first, stable whole-site roles for allele count tables."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re


REQUIRED={'cell','CHROM','POS','REF','ALT','AD','DP','GT','cM','gene'}


def canonical_site(chromosome,position):
    chromosome=str(chromosome)
    chromosome=chromosome if chromosome.startswith('chr') else 'chr'+chromosome
    if not re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)',chromosome):
        raise ValueError('Unsupported chromosome identity')
    if not re.fullmatch(r'[0-9]+',str(position)) or not 1<=int(position)<=2**31-1:
        raise ValueError('SNP POS needs a positive supported one-based integer')
    return chromosome,int(position)


def site_role(chromosome,position,*,genome_build,seed):
    chromosome,position=canonical_site(chromosome,position)
    if not genome_build or '\0' in genome_build:
        raise ValueError('Explicit genome-build identity required')
    payload=f'haplodose-allele-role-v1\0{int(seed)}\0{genome_build}\0{chromosome}\0{position}'.encode()
    value=int.from_bytes(hashlib.blake2b(payload).digest()[:8],'big')
    return 'discovery' if value<2**63 else 'reserved'


def split_table(path,out,*,genome_build,seed):
    """Partition before AD/DP, phase, coverage, bias or SNP-order filtering.

    Requires biallelic REF/ALT metadata in this research interface. Existing
    public readers still support their historical input format. Count strings
    are preserved; each partition is subsequently validated by the count reader.
    """
    path,out=Path(path),Path(out)
    out.mkdir(exist_ok=False)
    opener=gzip.open if path.suffix=='.gz' else open
    counts=dict(discovery=0,reserved=0);sites={};annotations={}
    with opener(path,'rt',newline='') as source, (out/'discovery.tsv').open('x',newline='') as discovery, (out/'reserved.tsv').open('x',newline='') as reserved:
        reader=csv.DictReader(source,delimiter='\t')
        fields=reader.fieldnames or []
        if len(fields)!=len(set(fields)) or not REQUIRED.issubset(fields):
            raise ValueError('Unique allele headers and explicit REF/ALT metadata required')
        writers={name:csv.DictWriter(stream,fieldnames=fields,delimiter='\t',lineterminator='\n')
                 for name,stream in (('discovery',discovery),('reserved',reserved))}
        for writer in writers.values():
            writer.writeheader()
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError('Malformed allele row')
            chrom,pos=canonical_site(row['CHROM'],row['POS'])
            ref,alt=row['REF'].upper(),row['ALT'].upper()
            if ref not in 'ACGT' or alt not in 'ACGT' or len(ref)!=1 or len(alt)!=1 or ref==alt:
                raise ValueError('Explicit distinct biallelic DNA bases required')
            key=f'{chrom}:{pos}'
            signature=(ref,alt,row['cM'],row['gene'])
            if key in annotations and annotations[key]!=signature:
                raise ValueError('Conflicting variant/map/gene annotation at '+key)
            annotations[key]=signature
            if key not in sites:
                sites[key]=site_role(chrom,pos,genome_build=genome_build,seed=seed)
            role=sites[key]
            writers[role].writerow(row);counts[role]+=1
    report=dict(schema='whole_site_roles_v1',genome_build=genome_build,seed=int(seed),
                counts=counts,sites=dict(sorted(sites.items())),
                required_assumption='REF/ALT and assembly are supplied correctly; not verified against FASTA')
    (out/'roles.json').write_text(json.dumps(report,indent=2))
    return report
