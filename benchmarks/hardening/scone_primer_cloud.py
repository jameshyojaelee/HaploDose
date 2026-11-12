"""Fetch the checksum-identified primary supplement through the current PMC service."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from urllib.parse import parse_qs, urlparse
from urllib.request import Request, urlopen

from haplodose.provenance import start_run, finish_run, file_record

BASE='https://pmc-oa-opendata.s3.amazonaws.com/'
EXPECTED_MD5='ea454b5fbfccc83e77de13f7ce8c5fe2'


def fetch(url,cap):
    assert url.startswith(BASE)
    with urlopen(Request(url,headers={'User-Agent':'HaploDose-primary-primer-source'}),timeout=60) as response:
        assert urlparse(response.url).hostname=='pmc-oa-opendata.s3.amazonaws.com'
        raw=response.read(cap+1)
    assert len(raw)<=cap
    return raw


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    binary=Path('/gpfs/commons/home/jameslee/micromamba/bin/pdftotext')
    run=start_run(args.out,dict(script=Path(__file__),design=args.task/'PRIMER_CLOUD_DESIGN.md',pdftotext=binary),
                  dict(metadata_url=BASE+'PMC9812385.1/PMC9812385.1.json',expected_MD5=EXPECTED_MD5),
                  dict(scope='Primary source identity and primer-table access, not executed-library provenance'),mode='scone_primer_cloud')
    try:
        metadata=fetch(BASE+'PMC9812385.1/PMC9812385.1.json',1024**2)
        (args.out/'article_metadata.json').write_bytes(metadata);record=json.loads(metadata)
        assert record['pmcid']=='PMC9812385' and record['version']==1
        assert record['doi']=='10.1126/sciadv.abp8901'
        assert record['is_manuscript'] is False and record['is_retracted'] is False
        matches=[value for value in record['media_urls'] if urlparse(value).path.endswith('/sciadv.abp8901_sm.pdf')]
        assert len(matches)==1;source=urlparse(matches[0])
        assert source.scheme=='s3' and source.netloc=='pmc-oa-opendata'
        assert source.path=='/PMC9812385.1/sciadv.abp8901_sm.pdf'
        assert parse_qs(source.query)=={'md5':[EXPECTED_MD5]}
        url=BASE+source.path.lstrip('/');pdf=fetch(url,32*1024**2)
        assert pdf.startswith(b'%PDF-') and hashlib.md5(pdf).hexdigest()==EXPECTED_MD5
        path=args.out/'sciadv.abp8901_sm.pdf';path.write_bytes(pdf)
        version=subprocess.run([str(binary),'-v'],capture_output=True,text=True,check=True)
        text=args.out/'supplement.txt'
        subprocess.run([str(binary),'-layout',str(path),str(text)],check=True,timeout=60)
        pages=text.read_text().split('\f')
        selected=[dict(page=i+1,text=page) for i,page in enumerate(pages)
                  if any(term in page for term in ('Table S2','TableS2','TCATG','AACTT','ATCGT'))]
        assert selected,'No declared primer-table or motif text extracted'
        (args.out/'primer_pages.json').write_text(json.dumps(selected,indent=2)+'\n')
        result=dict(source_accessible=True,url=url,declared_MD5=EXPECTED_MD5,checksum_verified=True,
                    pdf=file_record(path),selected_pages=[row['page'] for row in selected],
                    pdftotext=dict(**file_record(binary),version=(version.stdout+version.stderr).strip()),
                    exact_library_chemistry_verified=False,executed_author_workflow_verified=False,
                    biological_null_specimens_verified=0)
        (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
