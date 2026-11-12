"""Bounded primary archive metadata; never follow measurement-file links."""
import argparse,json
from pathlib import Path
from urllib.request import Request,urlopen
from haplodose.provenance import start_run,finish_run

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
 urls={'DNA_archive':'https://zenodo.org/api/records/13984673','GSM7016419':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7016419&targ=self&form=text&view=full','GSM7817780':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7817780&targ=self&form=text&view=full'}
 record=start_run(a.out,dict(script=Path(__file__)),dict(urls=urls,maximum_response_bytes=2097152),dict(scope='Correspondence feasibility only; no biological validation or per-cell data access'),mode='defnd_correspondence')
 try:
  reports=[]
  for name,url in urls.items():
   try:
    with urlopen(Request(url,headers={'User-Agent':'HaploDose-metadata-review'}),timeout=40) as h:data=h.read(2097153)
    if len(data)>2097152:raise ValueError('Metadata response exceeds cap')
    (a.out/(name+'.txt')).write_bytes(data)
    if name=='DNA_archive':
     obj=json.loads(data);meta=obj['metadata']
     item=dict(name=name,url=url,title=meta.get('title'),description=meta.get('description'),doi=meta.get('doi'),updated=obj.get('updated'),files=[dict(key=f['key'],size=f.get('size'),checksum=f.get('checksum'),url=f.get('links',{}).get('self')) for f in obj.get('files',[])])
    else:
     text=data.decode();lines=[line for line in text.splitlines() if line.startswith(('!Sample_title','!Sample_source_name','!Sample_characteristics','!Sample_extract_protocol','!Sample_data_processing','!Sample_relation','!Sample_supplementary_file','!Sample_platform_id'))]
     if not lines:raise ValueError('Expected GEO sample fields absent')
     item=dict(name=name,url=url,metadata_lines=lines)
    reports.append(item)
   except Exception as e:reports.append(dict(name=name,url=url,error=repr(e)))
  (a.out/'report.json').write_text(json.dumps(dict(records=reports,measurement_links_followed=False,cohort_admitted=False,biological_validation_n=0,limit='Archive labels alone do not establish same-cell correspondence, assembly or DNA detection sensitivity'),indent=2))
  finish_run(a.out,record,backend='cpu')
 except Exception as e:finish_run(a.out,record,error=e);raise
if __name__=='__main__':main()
