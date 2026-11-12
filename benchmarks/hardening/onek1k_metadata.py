"""Bounded archive metadata only; never retrieve cohort measurements or calls."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from haplodose.provenance import start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    studies = ('GSE196830', 'GSE196829', 'GSE196735')
    record = start_run(args.out, {'script': Path(__file__)},
        {'studies': studies, 'max_bytes_per_response': 2 * 1024 * 1024},
        {'scope': 'Feasibility only; no donor, cell, genotype or CNA outcomes'},
        mode='onek1k_metadata')
    reports = []
    try:
        for accession in studies:
            url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?' + urlencode(
                dict(acc=accession, targ='self', form='text', view='quick'))
            try:
                with urlopen(Request(url, headers={
                        'User-Agent': 'HaploDose-metadata-feasibility'}), timeout=30) as response:
                    raw = response.read(2 * 1024 * 1024 + 1)
                if len(raw) > 2 * 1024 * 1024:
                    raise ValueError('Archive metadata exceeds bounded request')
                text = raw.decode('utf-8')
                if '^SERIES = ' + accession not in text:
                    raise ValueError('Not the requested GEO SOFT series record')
                # A self-series response must contain no per-sample data tables.
                if '^SAMPLE' in text or '!sample_table_begin' in text.lower():
                    raise ValueError('Unexpected sample data in metadata response')
                path = args.out / (accession + '.soft.txt')
                path.write_bytes(raw)
                lines = text.splitlines()
                fields = {key: [line.split(' = ', 1)[1] for line in lines
                    if line.startswith('!Series_' + key + ' = ')] for key in
                    ('title', 'summary', 'overall_design', 'type', 'platform_id',
                     'sample_id', 'supplementary_file', 'relation', 'pubmed_id')}
                reports.append(dict(accession=accession, url=url, status='metadata_available',
                    bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest(), fields=fields,
                    sample_accessions_are_not_biological_n=True))
            except Exception as error:
                reports.append(dict(accession=accession, url=url,
                    status='metadata_unresolved', error=repr(error)))
        (args.out / 'summary.json').write_text(json.dumps(dict(reports=reports,
            measurements_downloaded=False, supplementary_links_followed=False,
            biological_null_specimens_verified=0,
            decision='Not admitted; correspondence, DNA sensitivity and independence unresolved'), indent=2))
        print(json.dumps([{k: v for k, v in report.items() if k != 'fields'}
            for report in reports]), flush=True)
        finish_run(args.out, record, backend='cpu')
    except Exception as error:
        finish_run(args.out, record, error=error)
        raise


if __name__ == '__main__':
    main()
