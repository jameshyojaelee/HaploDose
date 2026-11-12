"""Inspect a locked metadata-only pilot for RNA/array donor correspondence."""
import argparse
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from haplodose.provenance import start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    selection = json.loads((args.task/'selection.json').read_text())
    inputs = dict(script=Path(__file__), selection=args.task/'selection.json')
    for series, record in selection['series'].items():
        inputs[series] = Path(record['source'])
    run = start_run(args.out, inputs, dict(max_bytes=2097152, sample_records=6),
                    dict(scope='Metadata pilot only; no measurement or cohort admission'),
                    mode='onek1k_crosswalk')
    try:
        reports = []
        for series, record in selection['series'].items():
            assert len(record['lexical_pilot']) == 3
            for accession in record['lexical_pilot']:
                url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?' + urlencode(
                    dict(acc=accession, targ='self', form='text', view='quick'))
                report = dict(series=series, accession=accession, url=url)
                try:
                    with urlopen(Request(url, headers={'User-Agent':'HaploDose-metadata-pilot'}), timeout=45) as response:
                        raw = response.read(2097153)
                    assert len(raw) <= 2097152, 'Metadata size cap exceeded'
                    text = raw.decode('utf-8')
                    assert '^SAMPLE = ' + accession in text, 'Wrong sample record'
                    assert '!sample_table_begin' not in text.lower(), 'Measurement table in response'
                    (args.out/(accession+'.soft.txt')).write_bytes(raw)
                    fields = {}
                    for line in text.splitlines():
                        if line.startswith('!Sample_') and ' = ' in line:
                            key, value = line.split(' = ', 1)
                            fields.setdefault(key.removeprefix('!Sample_'), []).append(value)
                    report.update(status='metadata_available', bytes=len(raw), fields=fields)
                except Exception as error:
                    report.update(status='unresolved', error=repr(error))
                reports.append(report)
        (args.out/'report.json').write_text(json.dumps(dict(samples=reports,
            measurements_downloaded=False, supplementary_links_followed=False,
            biological_null_specimens_verified=0, cohort_admitted=False), indent=2))
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error); raise


if __name__ == '__main__':
    main()
