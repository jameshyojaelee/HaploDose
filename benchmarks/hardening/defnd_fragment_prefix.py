"""Bounded header/prefix inspection without RNA counts or cluster interpretation."""
import argparse,bz2,hashlib,json,re,zlib
from pathlib import Path
from urllib.request import Request,urlopen
from haplodose.provenance import start_run,finish_run
CAP=2*1024*1024;OUT=1024*1024

def prefix(url,kind,lines):
 dec=zlib.decompressobj(31) if kind=='gz' else bz2.BZ2Decompressor();raw=bytearray();data=bytearray()
 with urlopen(Request(url,headers={'Range':f'bytes=0-{CAP-1}','User-Agent':'HaploDose-bounded-correspondence'}),timeout=40) as h:
  headers={k:h.headers.get(k) for k in ['Content-Range','Content-Length','ETag','Last-Modified']};status=h.status
  while len(raw)<CAP and len(data)<OUT and sum(bool(v) and (kind=='gz' or not v.startswith(b'#')) for v in bytes(data).split(b'\n')[:-1])<lines:
   chunk=h.read(min(16384,CAP-len(raw)))
   if not chunk:break
   raw.extend(chunk);data.extend(dec.decompress(chunk,max_length=OUT-len(data)))
 kept=[];n=0
 for v in bytes(data).split(b'\n')[:-1]:
  kept.append(v)
  if v and (kind=='gz' or not v.startswith(b'#')):n+=1
  if n==lines:break
 if n<lines:raise ValueError('Bounded prefix did not contain requested complete data lines')
 selected=b'\n'.join(kept)+b'\n'
 return selected,dict(url=url,http_status=status,headers=headers,compressed_bytes_read=len(raw),compressed_prefix_sha256=hashlib.sha256(raw).hexdigest(),retained_bytes=len(selected),retained_sha256=hashlib.sha256(selected).hexdigest(),full_file_checksum_verified=False)

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
 root=a.task.parents[1];metadata=root/'dev/defnd_correspondence_20261108/metadata_v1/report.json';source=json.loads(metadata.read_text());archive=source['records'][0]
 rec=start_run(a.out,dict(script=Path(__file__),design=a.task/'PREFIX_DESIGN.md',metadata=metadata),dict(compressed_cap=CAP,decompressed_cap=OUT,dna_prefix_lines=64),dict(scope='Bounded admission only; not biological truth or performance'),mode='defnd_header_gate')
 try:
  reports=[]
  for item in source['records'][1:]:
   accession=item['name'];row=dict(accession=accession)
   try:
    url=next(v.split(' = ',1)[1] for v in item['metadata_lines'] if v.startswith('!Sample_supplementary_file_1')).replace('ftp://','https://')
    header,identity=prefix(url,'gz',1);(a.out/(accession+'_rna_header.txt')).write_bytes(header)
    fields=header.decode().rstrip('\n').split('\t')
    # A DNA alphabet token is only a candidate identity; no parsing rule is adopted.
    tokens=[re.findall(r'[ACGTN]{12,}',field) for field in fields]
    candidate=[v[0] for v in tokens if len(v)==1]
    row['rna']=dict(identity=identity,columns=len(fields),columns_with_one_dna_token=len(candidate),unique_candidate_barcodes=len(set(candidate)),ambiguous_token_columns=sum(len(v)>1 for v in tokens),parser_validated=False)
    dnaurl=next(f['url'] for f in archive['files'] if f['key'].startswith(accession+'_'))
    data,didentity=prefix(dnaurl,'bz2',64);(a.out/(accession+'_dna_prefix.tsv')).write_bytes(data)
    dna=[line.split('\t') for line in data.decode().splitlines() if line and not line.startswith('#')]
    row['dna']=dict(identity=didentity,data_rows=len(dna),field_counts=sorted({len(v) for v in dna}),comments=[line for line in data.decode().splitlines() if line.startswith('#')],chromosomes=sorted({v[0] for v in dna}),candidate_barcode_overlap_rows=sum(len(v)>=4 and v[3] in set(candidate) for v in dna),coordinate_integer_rows=sum(len(v)>=3 and v[1].isdigit() and v[2].isdigit() for v in dna),barcode_lengths=sorted({len(v[3]) for v in dna if len(v)>=4}),multiplicity_values=sorted({v[4] for v in dna if len(v)>=5}),positive_interval_rows=sum(len(v)>=3 and v[1].isdigit() and v[2].isdigit() and int(v[2])>int(v[1]) for v in dna))
   except Exception as e:row['unresolved']=repr(e)
   reports.append(row)
  (a.out/'report.json').write_text(json.dumps(dict(records=reports,cohort_admitted=False,RNA_count_rows_inspected=0,cluster_labels_used=False,biological_validation_n=0),indent=2))
  finish_run(a.out,rec,backend='cpu')
 except Exception as e:finish_run(a.out,rec,error=e);raise
if __name__=='__main__':main()
