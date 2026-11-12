"""Fetch only ENA read-run metadata for prospective feasibility, never reads/counts."""
import argparse
import csv
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from haplodose.provenance import start_run,finish_run

BASE='https://www.ebi.ac.uk/ena/portal/api/'
FIELDS=['study_accession','secondary_study_accession','study_alias','sample_accession',
        'sample_alias','experiment_accession','experiment_alias','run_accession',
        'experiment_title','library_strategy','library_source','library_layout',
        'instrument_model','read_count','base_count','fastq_ftp','fastq_md5','fastq_bytes']


def fetch(url):
    # Bounded metadata response; the whitelist below never requests a FASTQ,
    # count matrix, fragment file, clone assignment or repository source file.
    if not url.startswith(BASE):raise ValueError('Nonmetadata endpoint refused')
    with urlopen(Request(url,headers={'User-Agent':'HaploDose-metadata-feasibility'}),timeout=60) as r:
        content=r.read(10*1024*1024+1)
    if len(content)>10*1024*1024:raise ValueError('Metadata exceeds bounded request')
    return content


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    a.out.mkdir(exist_ok=False)
    rec=start_run(a.out,{'script':Path(__file__)},dict(fields=FIELDS,
        studies={'DEFND':'GSE224149','Hybrid_BAG':'PRJNA773107'}),
        dict(scope='Metadata eligibility only; no method execution or biological truth acceptance'),mode='prospective_metadata')
    reports=[]
    try:
        fields_url=BASE+'returnFields?'+urlencode({'result':'read_run'})
        raw=fetch(fields_url);(a.out/'available_fields.tsv').write_bytes(raw)
        fields=list(csv.DictReader(io.StringIO(raw.decode()),delimiter='\t'))
        available={r['columnId'] for r in fields}
        if not set(FIELDS)<=available:raise ValueError('Requested ENA metadata fields unavailable')
        for name,query in [('DEFND','study_alias="GSE224149"'),('Hybrid_BAG','study_accession="PRJNA773107"')]:
            url=BASE+'search?'+urlencode(dict(result='read_run',query=query,fields=','.join(FIELDS),format='tsv',limit=0))
            try:
                content=fetch(url);path=a.out/(name+'.tsv');path.write_bytes(content)
                rows=list(csv.DictReader(io.StringIO(content.decode()),delimiter='\t'))
                if not rows or not set(FIELDS)<=rows[0].keys():raise ValueError('No matching metadata or schema mismatch')
                if len({r['run_accession'] for r in rows})!=len(rows):raise ValueError('Duplicate run IDs')
                reports.append(dict(study=name,url=url,status='metadata_available',runs=len(rows),
                    sample_accessions=len({r['sample_accession'] for r in rows}),
                    library_strategies=dict(Counter(r['library_strategy'] for r in rows)),
                    library_sources=dict(Counter(r['library_source'] for r in rows)),
                    manifest_sha256=hashlib.sha256(content).hexdigest(),
                    biological_n='Not inferred from run/sample counts; independent donor identities require explicit curation',
                    truth_status='Not evaluated; DNA/RNA correspondence, callable resolution and DNA-only truth construction unresolved'))
            except Exception as error:
                reports.append(dict(study=name,url=url,status='metadata_unresolved',error=repr(error)))
        (a.out/'summary.json').write_text(json.dumps(dict(reports=reports,raw_data_downloaded=False,
            expression_or_clone_outputs_inspected=False,eligibility='No cohort accepted for confirmatory use'),indent=2))
        print(json.dumps(reports),flush=True)
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
