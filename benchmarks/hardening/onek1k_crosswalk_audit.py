"""Rebuild administrative joins directly from preserved primary metadata."""
import argparse
from collections import defaultdict, Counter
import gzip
import json
from pathlib import Path
import re

from haplodose.provenance import start_run, finish_run, file_record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args(); a.out.mkdir(exist_ok=False)
    source = a.task/'array_metadata_v1'
    manifest = json.loads((source/'run_manifest.json').read_text())
    assert manifest['status'] == 'complete'
    for group in ('input','output'):
        for record in manifest[group]['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
    inputs = dict(script=Path(__file__), manifest=source/'run_manifest.json',
                  report=source/'report.json', records=source/'records.jsonl',
                  series=Path(manifest['input']['files']['series']['path']),
                  selection=a.task/'selection.json')
    for path in (a.task/'donors_v1').glob('*.txt.gz'): inputs[path.name] = path
    for path in source.glob('*.soft.txt'): inputs[path.name] = path
    run = start_run(a.out,inputs,dict(expected_sample_records=1104),
                    dict(scope='Administrative identity; not biological or assay identity'),
                    mode='onek1k_crosswalk_audit')
    try:
        expected = sorted(set(re.findall(r'^!Series_sample_id = (GSM\d+)',inputs['series'].read_text(),re.M)))
        assert len(expected) == 1104
        records = [json.loads(s) for s in inputs['records'].read_text().splitlines()]
        assert [r['accession'] for r in records] == expected
        assert all(r['status']=='available' for r in records), 'Incomplete metadata cannot qualify the crosswalk'
        lookup = defaultdict(list); unknown = []; counts = []
        for accession in expected:
            lines = inputs[accession+'.soft.txt'].read_text().splitlines()
            assert '^SAMPLE = '+accession in lines
            assert not any(s.lower()=='!sample_table_begin' for s in lines)
            prefix = '!Sample_characteristics_ch1 = individual id: '
            identifiers = [s[len(prefix):] for s in lines if s.startswith(prefix)]
            if len(identifiers)==1 and identifiers[0]: lookup[identifiers[0]].append(accession)
            else: unknown.append(accession)
        selected = json.loads(inputs['selection'].read_text())['series']['GSE196735']['lexical_pilot']
        wanted = []
        for accession in selected:
            with gzip.open(inputs[accession+'.txt.gz'],'rt') as stream:
                values = stream.read().splitlines()
            assert all(re.fullmatch(r'[A-Za-z0-9_.-]{1,200}',s) for s in values)
            wanted.extend(values);counts.append(dict(accession=accession,identifiers=len(values)))
        assert len(wanted)==len(set(wanted))==40
        expected_matches = [dict(identifier=s, accessions=lookup[s], status=
            'unmatched' if not lookup[s] else 'unique' if len(lookup[s])==1 else 'ambiguous') for s in sorted(wanted)]
        reported = json.loads(inputs['report'].read_text())
        assert reported['attempted']==reported['expected']==1104 and reported['metadata_complete']
        assert reported['matches']==expected_matches
        assert reported['duplicate_ids']=={s:v for s,v in lookup.items() if len(v)>1}
        assert reported['biological_null_specimens_verified']==0 and not reported['cohort_admitted']
        (a.out/'report.json').write_text(json.dumps(dict(primary_metadata_records_checked=1104,
            donor_lists=counts,requested_ids_checked=40,join_status_counts=dict(Counter(s['status'] for s in expected_matches)),
            missing_or_multiple_individual_id_fields=unknown,administrative_duplicate_ids={s:v for s,v in lookup.items() if len(v)>1},
            biological_null_specimens_verified=0,cohort_admitted=False),indent=2))
        finish_run(a.out,run,backend='cpu')
    except Exception as error:
        finish_run(a.out,run,error=error);raise


if __name__ == '__main__':
    main()
