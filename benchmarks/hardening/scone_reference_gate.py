"""Read-only sequence, annotation, index and tool identity gate before raw scONE mapping."""
import argparse
from collections import Counter
import csv
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import struct
import subprocess

from haplodose.provenance import start_run, finish_run, file_record

REFERENCE=Path('/gpfs/commons/home/jameslee/CRC/data/external/ensembl112')
PARENT=Path('/gpfs/commons/home/slee/slee/reference_genome/GRCh38/refdata-gex-GRCh38-2020-A/fasta/genome.fa')
TOOLS=Path('/gpfs/commons/home/jameslee/micromamba/envs/rnaseq')
CANONICAL={str(i) for i in range(1,23)}|{'X','Y','MT'}


def rows(path):
    with path.open() as handle:return list(csv.DictReader(handle,delimiter='\t'))


def normalize(name):
    if name=='chrM':return 'MT'
    return name[3:] if name.startswith('chr') and name[3:] in CANONICAL else name


def static_elf(path):
    # A nonzero ldd status is legitimate only for a verified static ELF.
    with path.open('rb') as handle:
        header=handle.read(64);assert header[:4]==b'\x7fELF' and header[4] in (1,2) and header[5] in (1,2)
        endian='<' if header[5]==1 else '>'
        bits=32 if header[4]==1 else 64
        offset=struct.unpack_from(endian+('I' if bits==32 else 'Q'),header,28 if bits==32 else 32)[0]
        size,count=struct.unpack_from(endian+'HH',header,42 if bits==32 else 54)
        assert size>=4 and 0<count<1024
        types=[]
        for i in range(count):
            handle.seek(offset+i*size);types.append(struct.unpack(endian+'I',handle.read(4))[0])
    return dict(bits=bits,program_header_types=types,static=3 not in types and 2 not in types)


def fasta(path,parent=False):
    digest=hashlib.sha256();records={};name=None;length=0;sequence=hashlib.md5()
    def store():
        if name is not None:
            assert name not in records and length>0
            records[name]=dict(length=length,sequence_md5=sequence.hexdigest())
    with path.open('rb') as handle:
        for line in handle:
            digest.update(line)
            if line.startswith(b'>'):
                store();name=line[1:].split()[0].decode('ascii')
                if parent:name=normalize(name)
                length=0;sequence=hashlib.md5()
            else:
                bases=line.rstrip(b'\r\n');assert name is not None and not re.search(rb'[^ACGTNacgtn]',bases)
                length+=len(bases);sequence.update(bases.upper())
    store();assert CANONICAL<=set(records)
    return dict(path=str(path),bytes=path.stat().st_size,sha256=digest.hexdigest()),records


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    index=REFERENCE/'star_index';refname='Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    inputs=dict(script=Path(__file__),design=args.task/'REFERENCE_DESIGN.md',manifest=REFERENCE/'manifest.tsv',
                index_manifest=index/'index_files.sha256.tsv',build=index/'build.manifest.tsv',
                builder=Path('/gpfs/commons/home/jameslee/CRC/slurm/setup_reference.sbatch'))
    run=start_run(args.out,inputs,{},dict(scope='Reference and tool identity only; no alignments or biological validation'),mode='scone_reference_gate')
    try:
        with (REFERENCE.parent/'.setup_reference.lock').open('rb') as lock:
            fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
            metadata_before={p.name:file_record(p)['sha256'] for p in sorted((TOOLS/'conda-meta').glob('*')) if p.is_file()}
            manifest={row['file']:row for row in rows(inputs['manifest'])}
            ref_record,contigs=fasta(REFERENCE/refname)
            parent_record,parent_contigs=fasta(PARENT,parent=True)
            assert ref_record['sha256']==manifest[refname]['local_sha256']
            assert ref_record['bytes']==int(manifest[refname]['local_bytes'])
            assert parent_record['sha256']==manifest[refname]['remote_sha256']
            assert parent_record['bytes']==int(manifest[refname]['remote_bytes'])
            assert contigs==parent_contigs,'Normalization changed sequence or contig membership'
            verified={refname:ref_record}
            for name,row in manifest.items():
                if name==refname:continue
                record=file_record(REFERENCE/name)
                assert record['sha256']==row['local_sha256'] and record['bytes']==int(row['local_bytes'])
                verified[name]=record
            lengths={name:value['length'] for name,value in contigs.items()}
            fai={line.split('\t')[0]:int(line.split('\t')[1]) for line in (REFERENCE/(refname+'.fai')).read_text().splitlines()}
            assert lengths==fai
            star_lengths={line.split()[0]:int(line.split()[1]) for line in (index/'chrNameLength.txt').read_text().splitlines()}
            assert star_lengths==lengths
            build={row['key']:row['value'] for row in rows(inputs['build'])}
            gtfname='Homo_sapiens.GRCh38.112.gtf.gz'
            assert build['ready']=='TRUE' and build['reference_sha256']==ref_record['sha256']
            assert build['gtf_sha256']==verified[gtfname]['sha256']
            assert build['index_file_manifest_sha256']==file_record(inputs['index_manifest'])['sha256']
            index_records=[];listed=rows(inputs['index_manifest']);assert len(listed)==15
            for item in listed:
                path=Path(item['artifact']);assert path.parent.resolve()==index.resolve() and path.is_file()
                record=file_record(path);assert record['sha256']==item['sha256'] and record['bytes']==int(item['bytes'])
                index_records.append(record)
            annotations=Counter();absent=Counter();headers=[]
            with gzip.open(REFERENCE/gtfname,'rt') as handle:
                for line in handle:
                    if line.startswith('#'):
                        headers.append(line.rstrip());continue
                    fields=line.rstrip('\n').split('\t');assert len(fields)==9
                    chrom,kind=fields[0],fields[2];start,end=map(int,fields[3:5]);assert 1<=start<=end
                    if chrom not in lengths:absent[chrom]+=1;continue
                    assert end<=lengths[chrom],f'Annotation outside {chrom} reference length'
                    annotations[(chrom,kind)]+=1
            assert all(annotations[(chrom,'exon')]>0 for chrom in CANONICAL)
            (args.out/'reference_identity.json').write_text(json.dumps(dict(reference=verified,parent=parent_record,
                contigs=contigs,index_build=build,index_artifacts=index_records,
                annotation_records_on_absent_contigs=dict(absent),
                scope='Reference/annotation checks complete; tool gate and whole run not yet complete'),indent=2)+'\n')
            tools={}
            tool_env={**os.environ,'LD_LIBRARY_PATH':str(TOOLS/'lib')}
            for name,arguments,package in [('STAR',['--version'],'star'),('samtools',['--version'],'samtools'),
                                            ('featureCounts',['-v'],'subread')]:
                binary=TOOLS/'bin'/('STAR-plain' if name=='STAR' else name);record=file_record(binary)
                result=subprocess.run([str(binary),*arguments],env=tool_env,capture_output=True,text=True,timeout=30)
                assert result.returncode==0,f'{name} version failed: {result.stderr}'
                dependency=subprocess.run(['/usr/bin/ldd',str(binary)],env=tool_env,capture_output=True,text=True,timeout=30)
                elf=static_elf(binary)
                if dependency.returncode:
                    assert dependency.returncode==1 and elf['static'],'Dependency inspection failed for a dynamic executable'
                assert 'not found' not in dependency.stdout+dependency.stderr
                libraries={Path(match).resolve() for match in re.findall(r'(/[^\s()]+)',dependency.stdout)}
                package_paths=list((TOOLS/'conda-meta').glob(package+'-*.json'));assert len(package_paths)==1
                tools[name]=dict(binary=record,argv=[str(binary),*arguments],stdout=result.stdout,stderr=result.stderr,
                                 launcher=file_record(TOOLS/'bin'/name),
                                 LD_LIBRARY_PATH=tool_env['LD_LIBRARY_PATH'],ldd=dependency.stdout,
                                 ldd_stderr=dependency.stderr,ldd_exit_code=dependency.returncode,ELF=elf,
                                 libraries=[file_record(path) for path in sorted(libraries)],
                                 conda_record=json.loads(package_paths[0].read_text()))
            metadata_after={p.name:file_record(p)['sha256'] for p in sorted((TOOLS/'conda-meta').glob('*')) if p.is_file()}
            assert metadata_before==metadata_after,'Conda source environment changed during inventory'
            result=dict(reference=verified,parent=parent_record,contigs=contigs,
                        sequence_identity_after_primary_name_normalization=True,FAI_STAR_contigs_exact=True,
                        index_build=build,index_artifacts=index_records,annotation_headers=headers,
                        annotation_records=[dict(chromosome=c,feature=k,records=n) for (c,k),n in sorted(annotations.items())],
                        annotation_records_on_absent_contigs=dict(absent),tools=tools,conda_meta_inventory=metadata_after,
                        source_manifest_records=[file_record(path) for path in inputs.values()],
                        isolated_alignment_environment_ready=False,alignment_run=False,biological_null_specimens_verified=0)
            (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
            assert all(file_record(path)['sha256']==run['input']['files'][name]['sha256'] for name,path in inputs.items())
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
