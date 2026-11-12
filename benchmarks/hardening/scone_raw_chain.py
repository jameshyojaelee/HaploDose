"""One metadata-selected raw library: archive identity and bounded motif inspection."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import itertools
import json
from pathlib import Path
from urllib.parse import urlencode,urlparse
from urllib.request import Request,urlopen

from haplodose.provenance import start_run,finish_run,file_record

LIMIT=250000
MOTIFS={'DNA':'TCATG','RNA_author_FASTA':'AACTT','RNA_GEO_protocol':'ATCGT',
        'DNA_reverse_complement':'CATGA','RNA_FASTA_reverse_complement':'AAGTT',
        'RNA_GEO_reverse_complement':'ACGAT'}


def fastq(stream):
    while True:
        name=stream.readline()
        if not name:return
        sequence=stream.readline().rstrip(b'\r\n')
        plus=stream.readline();quality=stream.readline().rstrip(b'\r\n')
        assert name.startswith(b'@') and plus.startswith(b'+') and sequence
        assert len(sequence)==len(quality)
        token=name[1:].split()[0]
        if token.endswith((b'/1',b'/2')):token=token[:-2]
        yield token,sequence.upper()


def preflight():
    valid=b'@a/1\nTCATGAAAA\n+\nIIIIIIIII\n@b/1\nAACTTCCCC\n+\nIIIIIIIII\n'
    records=list(fastq(io.BytesIO(valid)))
    assert records==[(b'a',b'TCATGAAAA'),(b'b',b'AACTTCCCC')]
    invalid=[b'@a\nACTG\n+\nIII\n',b'@a\nACTG\n',b'a\nACTG\n+\nIIII\n']
    for raw in invalid:
        try:list(fastq(io.BytesIO(raw)))
        except AssertionError:pass
        else:raise AssertionError('Invalid FASTQ accepted')
    return dict(valid_records=2,invalid_FASTQ_controls=3)


def object_check(record):
    raw=record['content'].encode()
    assert hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==record['git_blob_sha']


def download(url,path,expected_bytes,expected_md5):
    assert url.startswith('https://ftp.sra.ebi.ac.uk/vol1/fastq/')
    md5=hashlib.md5();size=0
    with urlopen(Request(url,headers={'User-Agent':'HaploDose-primary-chain'}),timeout=45) as response,path.open('xb') as out:
        assert urlparse(response.url).hostname=='ftp.sra.ebi.ac.uk'
        while True:
            block=response.read(1024*1024)
            if not block:break
            size+=len(block);assert size<=expected_bytes
            md5.update(block);out.write(block)
    assert size==expected_bytes and md5.hexdigest()==expected_md5
    return dict(url=url,expected_bytes=expected_bytes,md5=md5.hexdigest(),**file_record(path))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True);p.add_argument('--metadata-task',type=Path,required=True)
    p.add_argument('--library-class',choices=['legacy_whole','later_nucleus'],default='legacy_whole')
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    design='DESIGN.md' if a.library_class=='legacy_whole' else 'NUCLEUS_DESIGN.md'
    inputs=dict(script=Path(__file__),design=a.task/design,barcodes=a.task/'barcodes.json',driver=a.task/'driver.json',
                author_sources=a.metadata_task/'source_records.json',ENA=a.metadata_task/'metadata_v1/ENA_runs.tsv',
                metadata_audit=a.metadata_task/'audit_v1/run_manifest.json')
    run=start_run(a.out,inputs,dict(library_class=a.library_class,limit_pairs=LIMIT,motifs=MOTIFS,max_compressed_bytes=2*1024**3),
                  dict(role='Development primary-chain feasibility, not prospective evaluation or biological validation'),mode='scone_raw_chain')
    try:
        controls=preflight()
        manifest=json.loads(inputs['metadata_audit'].read_text());assert manifest['status']=='complete'
        for record in manifest['output']['files'].values():
            assert file_record(record['path'])['sha256']==record['sha256']
        for path in (inputs['barcodes'],inputs['driver']):object_check(json.loads(path.read_text()))
        for record in json.loads(inputs['author_sources'].read_text())['records']:object_check(record)
        with inputs['ENA'].open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
        title_pattern='HCT116 scONEseq ' if a.library_class=='legacy_whole' else 'HCT116, extracted nuclei, cell '
        eligible=[r for r in rows if title_pattern in r['experiment_title'] and r['library_strategy']=='OTHER']
        assert eligible
        selected=min(eligible,key=lambda r:(r['sample_alias'],r['run_accession']))
        expected=('GSM5609521','SRR16194948') if a.library_class=='legacy_whole' else ('GSM6683025','SRR22047451')
        assert (selected['sample_alias'],selected['run_accession'])==expected
        assert selected['library_layout']=='PAIRED'
        urls=selected['fastq_ftp'].split(';');sizes=list(map(int,selected['fastq_bytes'].split(';')))
        hashes=selected['fastq_md5'].split(';');assert len(urls)==len(sizes)==len(hashes)==2
        assert sum(sizes)<=2*1024**3
        (a.out/'selected_library.json').write_text(json.dumps(dict(selection='Lexical GSM then run among OTHER titles containing '+title_pattern,
                    eligible_libraries=len(eligible),record=selected),indent=2)+'\n')
        url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?'+urlencode(dict(acc=selected['sample_alias'],targ='self',form='text',view='quick'))
        with urlopen(Request(url,headers={'User-Agent':'HaploDose-primary-chain'}),timeout=45) as response:
            raw=response.read(1024**2+1)
        text=raw.decode();assert len(raw)<=1024**2 and '^SAMPLE = '+selected['sample_alias'] in text
        assert '!sample_table_begin' not in text.lower()
        (a.out/'sample.soft.txt').write_bytes(raw)
        folder=a.out/'raw';folder.mkdir();files=[]
        for source,size,md5 in zip(urls,sizes,hashes):
            url='https://'+source
            name=Path(urlparse(url).path).name
            assert name.startswith(selected['run_accession']+'_') and name.endswith('.fastq.gz')
            files.append(download(url,folder/name,size,md5))
        prefix=[Counter(),Counter()];lengths=[Counter(),Counter()]
        offsets=[{name:Counter() for name in MOTIFS} for _ in range(2)]
        pairs=Counter();examined=0;headers=[]
        with gzip.open(files[0]['path'],'rb') as first,gzip.open(files[1]['path'],'rb') as second:
            iterator=itertools.islice(itertools.zip_longest(fastq(first),fastq(second)),LIMIT)
            for left,right in iterator:
                assert left is not None and right is not None and left[0]==right[0]
                if len(headers)<3:headers.append(left[0].decode())
                sequences=[left[1],right[1]];labels=[]
                for mate,seq in enumerate(sequences):
                    key=seq[:5].decode();prefix[mate][key]+=1;lengths[mate][len(seq)]+=1
                    labels.append(next((name for name,motif in MOTIFS.items() if key==motif),'other_prefix'))
                    for name,motif in MOTIFS.items():
                        position=seq.find(motif.encode(),0,35)
                        while position>=0:
                            offsets[mate][name][position]+=1
                            position=seq.find(motif.encode(),position+1,35)
                pairs['|'.join(labels)]+=1;examined+=1
        assert examined>0
        for mate in range(2):
            assert sum(prefix[mate].values())==sum(lengths[mate].values())==examined
            for name,motif in MOTIFS.items():assert offsets[mate][name][0]==prefix[mate][motif]
        assert sum(pairs.values())==examined
        report=dict(selected_sample=selected['sample_alias'],selected_run=selected['run_accession'],files=files,
                    pairs_examined=examined,maximum_pairs=LIMIT,first_three_pair_ids=headers,
                    prefix_5bp_counts=[dict(sorted(value.items())) for value in prefix],
                    read_length_counts=[dict(sorted(value.items())) for value in lengths],
                    motif_occurrences_start_0_through_30=[{name:dict(sorted(value.items())) for name,value in part.items()} for part in offsets],
                    paired_prefix_category_counts=dict(sorted(pairs.items())),controls=controls,
                    role='One development library, deterministic sequencing prefix; not random sampling or biological replication',
                    no_alignment=True,no_CNA_fit=True,no_published_labels=True,no_protocol_chosen=True,
                    warning='Exact motif counts are diagnostics, not classified molecules; observed prefixes cannot prove the archived executed workflow',
                    biological_null_specimens_verified=0)
        (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        (a.out/'output_files.json').write_text(json.dumps({str(f.relative_to(a.out)):file_record(f) for f in sorted(a.out.rglob('*'))
                    if f.is_file() and f.name not in ('run_manifest.json','output_files.json')},indent=2)+'\n')
        finish_run(a.out,run,backend='cpu')
    except Exception as error:
        finish_run(a.out,run,error=error);raise


if __name__=='__main__':main()
