"""Bounded identifier-list inspection; no RNA counts or DNA intensity retrieval."""
import argparse
from collections import Counter
import gzip
import io
import json
from pathlib import Path
import re
from urllib.request import Request, urlopen

from haplodose.provenance import start_run, finish_run


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args(); a.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), gate=a.task/'NEXT_GATE.md',
                  metadata=a.task/'metadata_v1/report.json')
    metadata = json.loads(inputs['metadata'].read_text())
    run = start_run(a.out, inputs, dict(compressed_cap=262144, decoded_cap=1048576),
                    dict(scope='Identifier metadata only; no cohort admission'), mode='onek1k_donor_lists')
    try:
        dna_ids = []
        for sample in metadata['samples']:
            assert sample['status'] == 'metadata_available'
            if sample['series'] == 'GSE196829':
                values = [s.split(': ',1)[1] for s in sample['fields']['characteristics_ch1']
                          if s.startswith('individual id: ')]
                assert len(values) == 1
                dna_ids.append(values[0])
        assert len(dna_ids) == len(set(dna_ids)) == 3
        reports = []; all_ids = []
        for sample in metadata['samples']:
            if sample['series'] != 'GSE196735': continue
            links = [v for k, values in sample['fields'].items() if k.startswith('supplementary_file')
                     for v in values if v.endswith('GenotypeSamples.txt.gz')]
            assert len(links) == 1
            original = links[0]
            assert original.startswith('ftp://ftp.ncbi.nlm.nih.gov/geo/')
            url = original.replace('ftp://','https://',1)
            report = dict(accession=sample['accession'], deposited_url=original, requested_url=url)
            try:
                with urlopen(Request(url, headers={'User-Agent':'HaploDose-identifier-metadata'}), timeout=45) as response:
                    raw = response.read(262145)
                    report.update(http_status=response.status, final_url=response.url)
                assert len(raw) <= 262144, 'Compressed cap exceeded'
                (a.out/(sample['accession']+'.txt.gz')).write_bytes(raw)
                with gzip.GzipFile(fileobj=io.BytesIO(raw)) as stream:
                    decoded = stream.read(1048577)
                assert len(decoded) <= 1048576, 'Decoded cap exceeded'
                text = decoded.decode('utf-8')
                lines = text.splitlines()
                # The allowed format is one literal ID per nonempty line.
                # Other formats remain unresolved without inspecting genotype values.
                if not lines or not all(re.fullmatch(r'[A-Za-z0-9_.-]{1,200}', s) for s in lines):
                    report.update(status='format_unresolved', nonempty_lines=sum(bool(s) for s in lines))
                    reports.append(report); continue
                counts = Counter(lines); all_ids.extend(lines)
                report.update(status='identifier_list', lines=len(lines), identifiers=lines,
                    duplicate_ids={k:v for k,v in counts.items() if v>1},
                    exact_pilot_array_matches=sorted(set(lines)&set(dna_ids)),
                    not_in_pilot_arrays=sorted(set(lines)-set(dna_ids)))
            except Exception as error:
                report.update(status='retrieval_unresolved', error=repr(error))
            reports.append(report)
        counts = Counter(all_ids)
        (a.out/'report.json').write_text(json.dumps(dict(lists=reports, pilot_array_ids=dna_ids,
            repeated_ids_across_lists={k:v for k,v in counts.items() if v>1},
            biological_null_specimens_verified=0, cohort_admitted=False,
            scope='Literal administrative identity only; no CNA status or independent biological n'),indent=2))
        finish_run(a.out,run,backend='cpu')
    except Exception as error:
        finish_run(a.out,run,error=error);raise


if __name__ == '__main__':
    main()
