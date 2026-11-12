"""Deposited HCT116 cell-label metadata and matrix-header identity gate."""
import argparse
from collections import Counter
import csv
import gzip
import io
import json
from pathlib import Path
import re
import time
from urllib.request import Request, urlopen

from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    root = args.task.parent
    metadata = root/'scone_feasibility_20261109/metadata_v1'
    raw = root/'scone_raw_chain_20261109/chain_v1'
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
                  prior=metadata/'run_manifest.json', run_table=metadata/'ENA_runs.tsv',
                  series=metadata/'GSE185269.soft.txt', selected=raw/'selected_library.json',
                  preserved_sample=raw/'sample.soft.txt')
    run = start_run(args.out, inputs, dict(selection='All literal HCT116 scONEseq titles'),
        dict(scope='Administrative sample labels only, no measurements or biological identity certification'),
        mode='scone_cell_identity')
    try:
        prior = json.loads(inputs['prior'].read_text()); assert prior['status'] == 'complete'
        for record in prior['output']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        with inputs['run_table'].open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        selected = []
        for row in rows:
            match = re.fullmatch(r'.*?: (GSM\d+): HCT116 scONEseq (\d+) Homo sapiens OTHER', row['experiment_title'])
            if match:
                assert match.group(1) == row['sample_alias']
                selected.append((row['sample_alias'], row['run_accession'], match.group(2)))
        selected.sort(); assert len(selected) == 98 and len({r[0] for r in selected}) == 98
        (args.out/'selected.json').write_text(json.dumps(selected, indent=2)+'\n')
        folder = args.out/'sample_metadata'; folder.mkdir()
        samples = []; requests = []
        for accession, run_id, title_number in selected:
            url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+accession+'&targ=self&form=text&view=quick'
            try:
                if accession == 'GSM5609521':
                    content = inputs['preserved_sample'].read_bytes()
                    source = 'Preserved locked pilot metadata'
                else:
                    with urlopen(Request(url, headers={'User-Agent':'HaploDose-identity-metadata'}), timeout=45) as response:
                        content = response.read(262145); source = response.url
                    assert len(content) <= 262144
                    time.sleep(0.25)
                text = content.decode('utf-8')
                assert '^SAMPLE = '+accession in text
                assert '!sample_table_begin' not in text.lower()
                path = folder/(accession+'.soft.txt'); path.write_bytes(content)
                fields = {}
                for line in text.splitlines():
                    if line.startswith('!') and ' = ' in line:
                        key, value = line[1:].split(' = ',1)
                        fields.setdefault(key, []).append(value)
                description_numbers = []
                for value in fields.get('Sample_description', []):
                    match = re.fullmatch(r'HCT116 single cell (\d+)', value)
                    if match: description_numbers.append(match.group(1))
                titles = fields.get('Sample_title', [])
                assert titles == ['HCT116 scONEseq '+title_number]
                samples.append(dict(GSM=accession, SRR=run_id, title=titles[0],
                    title_number_literal=title_number, description_numbers_literal=description_numbers,
                    descriptions=fields.get('Sample_description', []), source_names=fields.get('Sample_source_name_ch1', []),
                    relations=fields.get('Sample_relation', [])))
                requests.append(dict(GSM=accession, url=url, source=source, status='metadata_available', file=file_record(path)))
            except Exception as error:
                requests.append(dict(GSM=accession, url=url, status='unresolved', error=repr(error)))
            (args.out/'requests.json').write_text(json.dumps(requests, indent=2)+'\n')
        (args.out/'samples.json').write_text(json.dumps(samples, indent=2)+'\n')
        # Inspect precisely the first decompressed line. Gzip/HTTP buffering may
        # fetch additional compressed bytes; no later row is parsed or retained.
        links = [line.split(' = ',1)[1] for line in inputs['series'].read_text().splitlines()
                 if line.startswith('!Series_supplementary_file = ')]
        header_files = ['GSE185269_scONEseq_Cellline_RNA_gene_counts_matrix.csv.gz',
                        'GSE185269_scONEseq_Cellline_scWGS_normalized_counts_matrix.txt.gz']
        headers = []
        for name in header_files:
            matches = [link for link in links if link.rsplit('/',1)[-1] == name]; assert len(matches) == 1
            url = matches[0].replace('ftp://','https://',1)
            try:
                with urlopen(Request(url, headers={'User-Agent':'HaploDose-label-header-only'}), timeout=90) as response:
                    resolved = response.url
                    with gzip.GzipFile(fileobj=response) as decoded:
                        header = decoded.readline(1048577)
                assert header.endswith(b'\n') and len(header) <= 1048576, 'Truncated or oversized header'
                text = header.decode('utf-8-sig')
                delimiter = ',' if name.endswith('.csv.gz') else '\t'
                labels = next(csv.reader(io.StringIO(text), delimiter=delimiter))
                assert len(labels) > 1
                path = args.out/(name+'.header.txt'); path.write_bytes(header)
                headers.append(dict(deposited_url=matches[0], requested_url=url, resolved_url=resolved,
                    status='header_available', delimiter=delimiter, labels=labels,
                    duplicate_labels={k:n for k,n in Counter(labels).items() if n>1}, file=file_record(path)))
            except Exception as error:
                headers.append(dict(requested_url=url, status='unresolved', error=repr(error)))
            (args.out/'headers.json').write_text(json.dumps(headers, indent=2)+'\n')
        full = len(samples) == len(selected)
        one_description = all(len(row['description_numbers_literal']) == 1 for row in samples)
        offsets = Counter(int(row['description_numbers_literal'][0])-int(row['title_number_literal'])
            for row in samples if len(row['description_numbers_literal']) == 1)
        result = dict(selected_libraries=len(selected), metadata_records_available=len(samples),
            all_selected_metadata_available=full, every_available_record_has_one_description_number=one_description,
            description_title_offsets_descriptive=dict(offsets),
            unique_description_numbers=len({row['description_numbers_literal'][0] for row in samples
                if len(row['description_numbers_literal']) == 1}),
            headers_available=sum(row['status']=='header_available' for row in headers),
            correspondence_rule_adopted=False, expression_DNA_values_parsed=False,
            biological_identity_certified=False, cohort_admitted=False, biological_null_specimens_verified=0)
        (args.out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__': main()
