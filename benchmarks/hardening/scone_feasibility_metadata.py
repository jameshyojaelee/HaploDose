"""Bounded scONE archive metadata and source identity, without opening measurements."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode,urlparse
from urllib.request import Request,urlopen

from haplodose.provenance import start_run,finish_run,file_record

BASE='https://www.ebi.ac.uk/ena/portal/api/'
FIELDS=['study_accession','secondary_study_accession','sample_accession','sample_alias',
        'experiment_accession','experiment_alias','run_accession','experiment_title',
        'library_strategy','library_source','library_layout','read_count','base_count',
        'fastq_ftp','fastq_md5','fastq_bytes']


def fetch(url):
    assert url.startswith(BASE) or url.startswith('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?')
    with urlopen(Request(url,headers={'User-Agent':'HaploDose-metadata-feasibility'}),timeout=45) as response:
        assert urlparse(response.url).hostname in ('www.ebi.ac.uk','www.ncbi.nlm.nih.gov')
        content=response.read(10*1024*1024+1)
    assert len(content)<=10*1024*1024
    return content


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
    run=start_run(a.out,dict(script=Path(__file__),design=a.task/'DESIGN.md',sources=a.task/'source_records.json'),
                  dict(study='PRJNA768428',GEO='GSE185269',fields=FIELDS),
                  dict(scope='Metadata/source gate only; no measurements, clone labels or cohort acceptance'),mode='scone_feasibility_metadata')
    try:
        records=json.loads((a.task/'source_records.json').read_text())['records']
        for record in records:
            raw=record['content'].encode();header=f'blob {len(raw)}\0'.encode()
            assert hashlib.sha1(header+raw).hexdigest()==record['git_blob_sha'],record['path']
        requests=[]
        fields_url=BASE+'returnFields?'+urlencode(dict(result='read_run'))
        raw=fetch(fields_url);(a.out/'available_fields.tsv').write_bytes(raw)
        available={r['columnId'] for r in csv.DictReader(io.StringIO(raw.decode()),delimiter='\t')}
        assert set(FIELDS)<=available
        url=BASE+'search?'+urlencode(dict(result='read_run',query='study_accession="PRJNA768428"',
                                        fields=','.join(FIELDS),format='tsv',limit=0))
        raw=fetch(url);(a.out/'ENA_runs.tsv').write_bytes(raw)
        rows=list(csv.DictReader(io.StringIO(raw.decode()),delimiter='\t'))
        assert rows and set(FIELDS)<=set(rows[0])
        assert len({r['run_accession'] for r in rows})==len(rows)
        requests.append(dict(url=url,file='ENA_runs.tsv',status='metadata_available'))
        url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?'+urlencode(dict(acc='GSE185269',targ='self',form='text',view='quick'))
        try:
            raw=fetch(url);text=raw.decode()
            assert '^SERIES = GSE185269' in text
            assert '!sample_table_begin' not in text.lower() and '!series_matrix_table_begin' not in text.lower()
            (a.out/'GSE185269.soft.txt').write_bytes(raw)
            requests.append(dict(url=url,file='GSE185269.soft.txt',status='metadata_available'))
        except Exception as error:
            requests.append(dict(url=url,status='metadata_unresolved',error=repr(error)))
        from collections import Counter
        report=dict(source_blobs_verified=len(records),runs=len(rows),
                    samples=len({r['sample_accession'] for r in rows}),
                    library_strategies=dict(Counter(r['library_strategy'] for r in rows)),
                    library_sources=dict(Counter(r['library_source'] for r in rows)),requests=requests,
                    biological_n='Not inferred from cells, library aliases or runs',
                    biological_null_specimens_verified=0,cohort_admitted=False,
                    measurements_opened=False,raw_reads_opened=False,clone_labels_opened=False,
                    next_gate='Explicit patient/cell RNA-DNA correspondence, assay-compatible RNA counts and independent DNA truth with uncertainty')
        (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        (a.out/'output_files.json').write_text(json.dumps({str(f.relative_to(a.out)):file_record(f) for f in sorted(a.out.iterdir())
                 if f.is_file() and f.name not in ('run_manifest.json','output_files.json')},indent=2)+'\n')
        finish_run(a.out,run,backend='cpu')
    except Exception as error:
        finish_run(a.out,run,error=error);raise


if __name__=='__main__':main()
