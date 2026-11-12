"""Complete exact donor-ID lookup from bounded GEO self-sample metadata."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import re
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from haplodose.provenance import start_run, finish_run


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args(); a.out.mkdir(exist_ok=False)
    root = a.task.parents[1]
    inputs = dict(script=Path(__file__), gate=a.task/'ARRAY_METADATA_GATE.md',
                  series=root/'dev/onek1k_feasibility_20261108/metadata_v1/GSE196829.soft.txt',
                  donor_lists=a.task/'donors_v1/report.json')
    accessions = sorted(set(re.findall(r'^!Series_sample_id = (GSM\d+)', inputs['series'].read_text(), re.M)))
    assert len(accessions) == 1104
    donor_lists = json.loads(inputs['donor_lists'].read_text())
    assert all(r['status'] == 'identifier_list' for r in donor_lists['lists'])
    wanted = sorted({s for r in donor_lists['lists'] for s in r['identifiers']})
    assert len(wanted) == 40
    for accession in accessions:
        cached = a.task/'metadata_v1'/(accession+'.soft.txt')
        if cached.is_file(): inputs['cached:'+accession] = cached
    run = start_run(a.out,inputs,dict(samples=1104,max_bytes=2097152,requests_per_second=1),
                    dict(scope='Administrative metadata only; no measurement/CNA status'),
                    mode='onek1k_array_metadata')
    rows = []; lookup = defaultdict(list); consecutive = 0; last = 0.
    def save():
        matches = [dict(identifier=identifier, accessions=lookup[identifier],
                        status='unmatched' if not lookup[identifier] else
                        'unique' if len(lookup[identifier]) == 1 else 'ambiguous') for identifier in wanted]
        (a.out/'report.json').write_text(json.dumps(dict(attempted=len(rows),expected=1104,
            metadata_complete=len(rows)==1104 and all(r['status']=='available' for r in rows),
            matches=matches, duplicate_ids={k:v for k,v in lookup.items() if len(v)>1},
            biological_null_specimens_verified=0,cohort_admitted=False),indent=2))
    try:
        with (a.out/'records.jsonl').open('x') as journal:
            for accession in accessions:
                url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?' + urlencode(
                    dict(acc=accession,targ='self',form='text',view='quick'))
                row = dict(accession=accession,url=url)
                try:
                    if 'cached:'+accession in inputs:
                        raw = inputs['cached:'+accession].read_bytes(); row['source']='preserved_pilot'
                    else:
                        time.sleep(max(0.,1.-(time.monotonic()-last))); last=time.monotonic()
                        with urlopen(Request(url,headers={'User-Agent':'HaploDose-array-metadata'}),timeout=30) as response:
                            raw=response.read(2097153)
                        row['source']='archive_request'
                    assert len(raw)<=2097152, 'Metadata cap exceeded'
                    text=raw.decode('utf-8')
                    assert '^SAMPLE = '+accession in text, 'Wrong sample'
                    assert '!sample_table_begin' not in text.lower(), 'Measurement table in metadata'
                    (a.out/(accession+'.soft.txt')).write_bytes(raw)
                    fields=defaultdict(list)
                    for line in text.splitlines():
                        if line.startswith('!Sample_') and ' = ' in line:
                            key,value=line.split(' = ',1);fields[key.removeprefix('!Sample_')].append(value)
                    ids=[s.split(': ',1)[1] for s in fields['characteristics_ch1'] if s.startswith('individual id: ')]
                    row.update(status='available',fields=dict(fields),individual_ids=ids)
                    if len(ids)==1 and ids[0]: lookup[ids[0]].append(accession)
                    else: row['identity_status']='missing_or_multiple_individual_ids'
                    consecutive=0
                except Exception as error:
                    row.update(status='unresolved',error=repr(error));consecutive+=1
                rows.append(row);journal.write(json.dumps(row)+'\n');journal.flush()
                if len(rows)%25==0:
                    save();print(json.dumps(dict(attempted=len(rows),total=1104)),flush=True)
                if consecutive>=3:
                    raise RuntimeError('Three consecutive metadata failures; partial records preserved')
        save();finish_run(a.out,run,backend='cpu')
    except Exception as error:
        save();finish_run(a.out,run,error=error);raise


if __name__ == '__main__':
    main()
