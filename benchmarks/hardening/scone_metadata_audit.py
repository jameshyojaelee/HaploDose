"""Independently reconcile scONE archive identifiers and record unresolved raw access."""
import argparse
import csv
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

from haplodose.provenance import start_run,finish_run,file_record


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
    data=a.task/'metadata_v1'
    run=start_run(a.out,dict(script=Path(__file__),design=a.task/'AUDIT_DESIGN.md',
                            metadata=data/'run_manifest.json',sources=a.task/'source_records.json'),{},
                  dict(scope='Administrative metadata, not specimen or true-clone validation'),mode='scone_metadata_audit')
    try:
        manifest=json.loads((data/'run_manifest.json').read_text());assert manifest['status']=='complete'
        for kind in ('input','output'):
            for record in manifest[kind]['files'].values():
                assert file_record(record['path'])['sha256']==record['sha256']
        with (data/'ENA_runs.tsv').open() as handle:
            rows=list(csv.DictReader(handle,delimiter='\t'))
        text=(data/'GSE185269.soft.txt').read_text()
        ids=re.findall(r'^!Series_sample_id = (GSM\d+)$',text,re.M)
        assert len(ids)==len(set(ids)) and all(r['sample_alias'] for r in rows)
        aliases={r['sample_alias'] for r in rows};geo=set(ids)
        report=json.loads((data/'report.json').read_text())
        assert len(rows)==report['runs'] and len({r['sample_accession'] for r in rows})==report['samples']
        assert dict(Counter(r['library_strategy'] for r in rows))==report['library_strategies']
        source_checks=[]
        for record in json.loads((a.task/'source_records.json').read_text())['records']:
            raw=record['content'].encode()
            # Independent binary object reconstruction, including retained CRLF.
            git=hashlib.sha1(b'blob '+str(len(raw)).encode()+bytes([0])+raw).hexdigest()
            assert git==record['git_blob_sha']
            source_checks.append(dict(path=record['path'],git_blob_sha=git,sha256=hashlib.sha256(raw).hexdigest()))
        titles=Counter();brain=[]
        for row in rows:
            title=row['experiment_title']
            match=re.search(r'GSM\d+: (.+?) Homo sapiens',title)
            name=match.group(1) if match else title
            titles[re.sub(r'\b\d+\b','#',name)]+=1
            if re.search(r'GBM|glioma|astro|brain|patient',name,re.I):
                brain.append(dict(run=row['run_accession'],sample_alias=row['sample_alias'],title=title))
        with (a.out/'title_patterns.tsv').open('w') as handle:
            writer=csv.writer(handle,delimiter='\t');writer.writerow(['literal_title_with_standalone_numbers_masked','runs'])
            writer.writerows(sorted(titles.items()))
        links=re.findall(r'^!Series_supplementary_file = (.+)$',text,re.M)
        (a.out/'supplementary_links.txt').write_text('\n'.join(links)+'\n')
        result=dict(runs=len(rows),GEO_sample_ids=len(ids),unique_ENA_sample_aliases=len(aliases),
                    GEO_without_ENA_alias=sorted(geo-aliases),ENA_alias_without_GEO=sorted(aliases-geo),
                    source_identity_checks=source_checks,literal_brain_title_matches=brain,
                    literal_brain_title_match_count=len(brain),title_patterns=len(titles),
                    deposited_GBM_links=[link for link in links if 'GBM' in link],
                    raw_brain_access='Not established by these records; a zero title-match count alone cannot establish absence',
                    metadata_reconciliation_only=True,cohort_admitted=False,biological_null_specimens_verified=0,
                    interpretation='No donor replication or DNA-RNA specimen identity inferred from equal administrative IDs')
        (a.out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
        finish_run(a.out,run,backend='cpu')
    except Exception as error:
        finish_run(a.out,run,error=error);raise


if __name__=='__main__':main()
