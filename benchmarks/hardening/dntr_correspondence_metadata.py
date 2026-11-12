"""Primary library metadata only: exact DNTR-seq RNA/DNA correspondence gate."""
import argparse
from collections import Counter,defaultdict
import csv
from decimal import Decimal,InvalidOperation
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
from urllib.parse import urlencode,urlparse
from urllib.request import Request,urlopen

from haplodose.provenance import start_run,finish_run,file_record

ENA='https://www.ebi.ac.uk/ena/portal/api/'
GEO='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144296/soft/GSE144296_family.soft.gz'
FIELDS=['run_accession','study_accession','sample_accession','experiment_accession','experiment_title',
        'library_strategy','library_source','library_layout','read_count','base_count','fastq_ftp','fastq_md5','fastq_bytes']
KEEP={'title','geo_accession','source_name_ch1','molecule_ch1','library_strategy','library_source',
      'library_selection','relation','data_processing','supplementary_file','taxid_ch1'}
CHARACTERISTICS={'cell id','cell line','cell type','cell count','treatment','time point'}
MISSING={'','na','n/a','unknown','not available','not provided','not applicable'}


def fetch(url):
    if not (url.startswith(ENA) or url==GEO):raise ValueError('Only declared primary metadata endpoints allowed')
    with urlopen(Request(url,headers={'User-Agent':'HaploDose-metadata-correspondence'}),timeout=45) as response:
        if urlparse(response.url).hostname not in ('www.ebi.ac.uk','ftp.ncbi.nlm.nih.gov'):raise ValueError('Unexpected metadata redirect')
        raw=response.read(10*1024*1024+1)
    if len(raw)>10*1024*1024:raise ValueError('Metadata exceeds compressed/response bound')
    return raw


def parse_geo(lines):
    records=[];current=None;kept=[];characters=Counter()
    for line in lines:
        if line.lower().startswith(('!sample_table_begin','!series_matrix_table_begin')):
            raise ValueError('Measurement table refused in metadata-only gate')
        if line.startswith('^'):
            current=None
            if line.startswith('^SAMPLE = '):
                accession=line.split(' = ',1)[1].strip();current=dict(accession=accession,fields={},characteristics={})
                records.append(current);kept.append(line.rstrip())
        elif current is not None and line.startswith('!Sample_') and ' = ' in line:
            key,value=line.rstrip().split(' = ',1);key=key.removeprefix('!Sample_')
            if key=='characteristics_ch1' and ': ' in value:
                name,content=value.split(': ',1);name=name.lower().strip();characters[name]+=1
                if name in CHARACTERISTICS:
                    if name in current['characteristics']:raise ValueError('Duplicate specimen characteristic')
                    current['characteristics'][name]=content;kept.append(line.rstrip())
            elif key in KEEP:
                current['fields'].setdefault(key,[]).append(value);kept.append(line.rstrip())
    if not records or len({r['accession'] for r in records})!=len(records):raise ValueError('Empty/duplicate GEO samples')
    return records,'\n'.join(kept)+'\n',dict(characters)


def role(record):
    molecules=record['fields'].get('molecule_ch1',[])
    if len(molecules)!=1:return 'unresolved'
    value=molecules[0].lower()
    if value=='genomic dna':return 'DNA'
    if value in ('total rna','poly a rna','polya rna'):return 'RNA'
    return 'unresolved'


def pairs(records):
    cells=defaultdict(list);unkeyed=[]
    for record in records:
        cell=record['characteristics'].get('cell id')
        if cell:cells[cell].append(record)
        else:unkeyed.append(record['accession'])
    output=[]
    for cell,rows in sorted(cells.items()):
        groups={label:[r for r in rows if role(r)==label] for label in ('RNA','DNA','unresolved')}
        status='paired_metadata_only';reasons=[]
        if len(groups['RNA'])!=1 or len(groups['DNA'])!=1 or groups['unresolved']:
            status='unresolved_library_cardinality';reasons.append({k:len(v) for k,v in groups.items()})
        else:
            for name in ('cell line','cell type','cell count','treatment','time point'):
                values=[r['characteristics'].get(name) for r in rows]
                if any(v is None or str(v).strip().casefold() in MISSING for v in values) or len(set(values))!=1:
                    reasons.append(dict(characteristic=name,values=values));status='discordant_or_missing_metadata'
            counts=[r['characteristics'].get('cell count','') for r in rows]
            try:single=all(Decimal(v)==1 for v in counts)
            except InvalidOperation:single=False
            if not single:status='not_certified_single_cell';reasons.append(dict(cell_counts=counts))
        output.append(dict(cell_id=cell,status=status,reasons=reasons,
            rna_samples=[r['accession'] for r in groups['RNA']],dna_samples=[r['accession'] for r in groups['DNA']],
            libraries=[dict(accession=r['accession'],role=role(r),characteristics=r['characteristics']) for r in rows]))
    return output,unkeyed


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--metadata-source',type=Path)
    a=p.parse_args();a.out.mkdir(exist_ok=False);root=a.task.parents[1]
    frozen=root/'benchmarks/heldout2/dntrseq_gse144296_runs.tsv';local=root/'data/external/dntrseq_gse144296'
    rec=start_run(a.out,dict(script=Path(__file__),design=a.task/'DESIGN.md',frozen_runs=frozen,
        download_manifest=local/'manifest.tsv',exposure_history=local/'README.txt'),dict(study='PRJNA603321',GEO='GSE144296',fields=FIELDS),
        dict(scope='Existing exposed regression cohort; administrative correspondence only,not independent validation or DNA truth'),mode='dntr_correspondence_metadata')
    requests=[]
    try:
        url=ENA+'search?'+urlencode(dict(result='read_run',query='study_accession="PRJNA603321"',fields=','.join(FIELDS),format='tsv',limit=0))
        if a.metadata_source:
            manifest=json.loads((a.metadata_source/'run_manifest.json').read_text());assert manifest['status']=='complete'
            for name in ('ENA_runs.tsv','GEO_metadata_only.soft'):
                record=manifest['output']['files'][name];assert file_record(record['path'])['sha256']==record['sha256']
            raw=(a.metadata_source/'ENA_runs.tsv').read_bytes()
            rec['input']['files']['cached_metadata_manifest']=file_record(a.metadata_source/'run_manifest.json')
            rec['input']['files']['cached_ENA']=file_record(a.metadata_source/'ENA_runs.tsv')
            rec['input']['files']['cached_GEO']=file_record(a.metadata_source/'GEO_metadata_only.soft')
        else:raw=fetch(url)
        (a.out/'ENA_runs.tsv').write_bytes(raw)
        ena=list(csv.DictReader(io.StringIO(raw.decode()),delimiter='\t'));assert ena and set(FIELDS)<=ena[0].keys()
        assert len({r['run_accession'] for r in ena})==len(ena)
        requests.append(dict(url=url,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
        if a.metadata_source:
            decoded=(a.metadata_source/'GEO_metadata_only.soft').read_bytes()
            requests.append(dict(url=GEO,reused_primary_metadata=str(a.metadata_source),sha256=hashlib.sha256(decoded).hexdigest()))
        else:
            raw=fetch(GEO);requests.append(dict(url=GEO,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
            decoded=gzip.GzipFile(fileobj=io.BytesIO(raw)).read(100*1024*1024+1)
        if len(decoded)>100*1024*1024:raise ValueError('Decompressed metadata exceeds bound')
        records,extracted,keys=parse_geo(io.StringIO(decoded.decode()))
        (a.out/'GEO_metadata_only.soft').write_text(extracted)
        (a.out/'GEO_samples.json').write_text(json.dumps(records,indent=2)+'\n')
        by_exp=defaultdict(list)
        for record in records:
            experiments=set(re.findall(r'\bSRX\d+\b',' '.join(record['fields'].get('relation',[]))))
            for experiment in experiments:by_exp[experiment].append(record)
        old=list(csv.DictReader(frozen.open(),delimiter='\t'));frozen_by_run={r['run']:r for r in old};assert len(frozen_by_run)==len(old)
        archived=list(csv.DictReader((local/'manifest.tsv').open(),delimiter='\t'));by_file={r['file']:r for r in archived};assert len(by_file)==len(archived)
        old_set=set(frozen_by_run);new_set={r['run_accession'] for r in ena};selected=[r for r in ena if r['run_accession'] in old_set]
        run_rows=[];file_rows=[];seen=set()
        for run in selected:
            previous=frozen_by_run[run['run_accession']];samples=by_exp[run['experiment_accession']]
            seen.update(r['accession'] for r in samples)
            run_rows.append(dict(run=run['run_accession'],experiment=run['experiment_accession'],samples=[r['accession'] for r in samples],
                explicit_cell_ids=[r['characteristics'].get('cell id') for r in samples],roles=[role(r) for r in samples],
                frozen_group=previous['group'],geo_join='unique' if len(samples)==1 else 'missing_or_ambiguous',
                read_metadata_unchanged=all(run[k]==previous[k] for k in ('read_count','base_count','fastq_md5','fastq_bytes'))))
            ftp,md5,size=[run[k].split(';') for k in ('fastq_ftp','fastq_md5','fastq_bytes')]
            assert len(ftp)==len(md5)==len(size)
            for address,checksum,n in zip(ftp,md5,size):
                relative=Path('fastq')/previous['group']/run['run_accession']/address.rsplit('/',1)[-1]
                record=by_file.get(relative.as_posix());path=local/relative
                archived_match=record is not None and record['md5']==checksum and record['bytes']==n
                exists=path.is_file();size_match=exists and path.stat().st_size==int(n)
                file_rows.append(dict(run=run['run_accession'],file=relative.as_posix(),archived_metadata_matches=archived_match,
                    exists=exists,size_matches=size_match,current_content_md5_recomputed=False))
        cell_ids={r['characteristics'].get('cell id') for r in records if r['accession'] in seen}
        cell_ids.discard(None)
        relevant=[r for r in records if r['accession'] in seen or r['characteristics'].get('cell id') in cell_ids]
        correspondence,unkeyed=pairs(relevant)
        for pair in correspondence:
            pair['outside_frozen_library_scope']=[x['accession'] for x in pair['libraries'] if x['accession'] not in seen]
            if pair['outside_frozen_library_scope'] and pair['status']=='paired_metadata_only':pair['status']='outside_frozen_library_scope'
        (a.out/'run_correspondence.json').write_text(json.dumps(run_rows,indent=2)+'\n')
        (a.out/'file_presence.json').write_text(json.dumps(file_rows,indent=2)+'\n')
        (a.out/'cell_correspondence.json').write_text(json.dumps(correspondence,indent=2)+'\n')
        report=dict(ena_runs=len(ena),frozen_runs=len(old),geo_samples=len(records),linked_geo_samples=len(relevant),
            missing_frozen_runs=sorted(old_set-new_set),additional_current_runs=sorted(new_set-old_set),
            run_join_status=dict(Counter(r['geo_join'] for r in run_rows)),cell_status=dict(Counter(r['status'] for r in correspondence)),
            samples_without_explicit_cell_id=unkeyed,read_metadata_changed=[r['run'] for r in run_rows if not r['read_metadata_unchanged']],
            expected_files=len(file_rows),missing_or_wrong_size=sum(not r['size_matches'] for r in file_rows),
            archive_metadata_disagreements=sum(not r['archived_metadata_matches'] for r in file_rows),
            characteristic_field_names=keys,current_content_checksums_recomputed=False,
            biological_specimens_added=0,certified_single_CNA_null_specimens=0,cohort_admitted=False,
            measurements_or_clone_labels_opened=False,raw_reads_opened=False,
            next_gate='Resolve explicit metadata conflicts; verify molecular correspondence and DNA callable sensitivity independently before regression scoring')
        (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        (a.out/'requests.json').write_text(json.dumps(requests,indent=2)+'\n')
        for record in rec['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(a.out,rec,backend='CPU bounded primary metadata and file-stat audit;no FASTQ/count reading')
    except BaseException as error:
        (a.out/'requests.json').write_text(json.dumps(dict(completed=requests,error=repr(error)),indent=2)+'\n')
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
