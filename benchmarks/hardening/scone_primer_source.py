"""Retrieve the primary primer supplement once; record unavailable access explicitly."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET
import zipfile

from haplodose.provenance import start_run, finish_run, file_record

URL='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9812385/supplementaryFiles'
NCBI_URL='https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC9812385'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--source',choices=['europepmc','ncbi'],default='europepmc')
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    ncbi=args.source=='ncbi'
    url=NCBI_URL if ncbi else URL
    cap=(64 if ncbi else 40)*1024**2
    design='PRIMER_NCBI_DESIGN.md' if ncbi else 'PRIMER_SOURCE_DESIGN.md'
    binary=Path('/gpfs/commons/home/jameslee/micromamba/bin/pdftotext')
    run=start_run(args.out,dict(script=Path(__file__),design=args.task/design,
                               pdftotext=binary),dict(url=url,maximum_bytes=cap),
                  dict(scope='Primary source access only; no count inference or biological validation'),mode='scone_primer_source')
    try:
        result=dict(url=url,source_accessible=False,primer_layout_verified=False,
                    biological_null_specimens_verified=0)
        try:
            if ncbi:
                with urlopen(Request(url,headers={'User-Agent':'HaploDose-primary-primer-source'}),timeout=60) as response:
                    metadata=response.read(1024**2+1)
                assert len(metadata)<=1024**2
                (args.out/'oa_record.xml').write_bytes(metadata)
                tree=ET.fromstring(metadata)
                records=tree.findall('.//record');assert len(records)==1 and records[0].get('id')=='PMC9812385'
                links=[link for link in records[0].findall('link') if link.get('format')=='tgz']
                assert len(links)==1
                url=links[0].get('href','').replace('ftp://','https://',1)
                assert url.startswith('https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_package/')
                result['package_url']=url
            with urlopen(Request(url,headers={'User-Agent':'HaploDose-primary-primer-source'}),timeout=60) as response:
                result.update(final_url=response.url,status=response.status,content_type=response.headers.get('Content-Type'))
                raw=response.read(cap+1)
        except (HTTPError,URLError,TimeoutError) as error:
            result['access_error']=str(error)
        else:
            assert len(raw)<=cap,'Supplement exceeds prespecified transport cap'
            (args.out/'response.bin').write_bytes(raw)
            result.update(response_bytes=len(raw),response_sha256=hashlib.sha256(raw).hexdigest())
            if not ncbi and not zipfile.is_zipfile(io.BytesIO(raw)):
                result['access_error']='Response is not the required supplementary ZIP'
            else:
                container=tarfile.open(fileobj=io.BytesIO(raw),mode='r:gz') if ncbi else zipfile.ZipFile(io.BytesIO(raw))
                with container as archive:
                    names=archive.getnames() if ncbi else archive.namelist()
                    assert len(names)==len(set(names))
                    result['archive_members']=names
                    matches=[name for name in names if Path(name).name=='sciadv.abp8901_sm.pdf']
                    if len(matches)!=1:
                        result['access_error']='Expected unique primary supplement PDF not found'
                    else:
                        if ncbi:
                            member=archive.getmember(matches[0]);assert member.isfile() and member.size<=32*1024**2
                            with archive.extractfile(member) as handle:pdf=handle.read(32*1024**2+1)
                        else:
                            member=archive.getinfo(matches[0]);assert member.file_size<=32*1024**2
                            pdf=archive.read(member)
                        assert len(pdf)<=32*1024**2 and pdf.startswith(b'%PDF-')
                        destination=args.out/'sciadv.abp8901_sm.pdf';destination.write_bytes(pdf)
                        version=subprocess.run([str(binary),'-v'],capture_output=True,text=True,check=True)
                        result['pdftotext']=dict(**file_record(binary),version=(version.stdout+version.stderr).strip())
                        text=args.out/'supplement.txt'
                        subprocess.run([str(binary),'-layout',str(destination),str(text)],check=True,timeout=60)
                        pages=text.read_text().split('\f')
                        selected=[dict(page=i+1,text=page) for i,page in enumerate(pages)
                                  if any(term in page for term in ('Table S2','TableS2','TCATG','AACTT','ATCGT'))]
                        (args.out/'primer_pages.json').write_text(json.dumps(selected,indent=2)+'\n')
                        result.update(source_accessible=True,pdf_member=matches[0],selected_pages=[row['page'] for row in selected])
        (args.out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
