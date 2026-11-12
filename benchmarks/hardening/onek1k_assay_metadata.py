"""Summarize assay-resource provenance from audited primary GEO metadata only."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import re

from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    source = args.task / 'array_metadata_v1'
    inputs = dict(script=Path(__file__), audit=args.task/'audit_v1/run_manifest.json',
                  report=source/'report.json', manifest=source/'run_manifest.json',
                  gate=args.task/'ASSAY_METADATA_GATE.md')
    manifests = [json.loads(inputs[key].read_text()) for key in ('audit', 'manifest')]
    for manifest in manifests:
        assert manifest['status'] == 'complete'
        for group in ('input', 'output'):
            for record in manifest[group]['files'].values():
                assert file_record(record['path'])['sha256'] == record['sha256']
    for path in source.glob('*.soft.txt'):
        inputs[path.name] = path
    run = start_run(args.out, inputs, {},
                    dict(scope='Metadata statements, not verified assay identity or CNA truth'),
                    mode='onek1k_assay_metadata')
    try:
        report = json.loads(inputs['report'].read_text())
        assert report['metadata_complete'] and report['attempted'] == 1104
        matched = {s for row in report['matches'] for s in row['accessions']}
        assert len(matched) == 40 and all(row['status'] == 'unique' for row in report['matches'])
        processing = Counter()
        matched_processing = Counter()
        platforms = Counter()
        pair_issues = []
        assays = defaultdict(list)
        rows = []
        for path in sorted(source.glob('*.soft.txt')):
            fields = defaultdict(list)
            for line in path.read_text().splitlines():
                if line.startswith('!Sample_') and ' = ' in line:
                    key, value = line.split(' = ', 1)
                    fields[key.removeprefix('!Sample_')].append(value)
            statements = fields['data_processing']
            versions = sorted({v for s in statements for v in re.findall(
                r'Global Screening Array-24 (v[0-9.]+)', s)})
            processing.update(statements)
            if path.stem.removesuffix('.soft') in matched:
                matched_processing.update(statements)
            platforms.update(fields['platform_id'])
            accession = path.name.removesuffix('.soft.txt')
            pairs = []
            for url in fields['supplementary_file']:
                match = re.fullmatch(r'(GSM\d+)_(\d+)_(R\d+C\d+)_(Grn|Red)\.idat\.gz', url.rsplit('/', 1)[-1])
                if match:
                    pairs.append(match.groups())
            valid = (len(pairs) == 2 and {p[3] for p in pairs} == {'Grn', 'Red'}
                     and len({p[:3] for p in pairs}) == 1
                     and all(p[0] == accession for p in pairs))
            if not valid:
                pair_issues.append(dict(accession=accession, links=fields['supplementary_file']))
            else:
                assays[':'.join(pairs[0][1:3])].append(accession)
            rows.append(dict(accession=accession, pilot_donor_match=accession in matched,
                             stated_versions=versions, pair_metadata_consistent=valid))
        assert len(rows) == 1104
        summary = dict(sample_records=len(rows), matched_pilot_array_records=len(matched),
                       data_processing_statements=dict(processing),
                       matched_data_processing_statements=dict(matched_processing),
                       version_statement_counts=dict(Counter(','.join(r['stated_versions']) or 'missing' for r in rows)),
                       platforms=dict(platforms), pair_metadata_issues=pair_issues,
                       repeated_chip_positions={k:v for k,v in assays.items() if len(v)>1},
                       assay_version_verified=False, conversion_resources_selected=False,
                       biological_null_specimens_verified=0, cohort_admitted=False)
        (args.out/'samples.json').write_text(json.dumps(rows, indent=2))
        (args.out/'report.json').write_text(json.dumps(summary, indent=2))
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
