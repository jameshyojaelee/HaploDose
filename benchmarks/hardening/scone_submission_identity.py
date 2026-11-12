"""Bounded original-submission metadata check for the locked scONE cell identities."""
import argparse
import csv
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    original = args.task.parent/'scone_feasibility_20261109/metadata_v1'
    inputs = dict(script=Path(__file__), design=args.task/'SUBMISSION_DESIGN.md',
        identity=args.task/'identity_v1/run_manifest.json',
        fields=original/'available_fields.tsv', earlier_runs=original/'ENA_runs.tsv')
    run = start_run(args.out, inputs, {}, dict(scope='Submitted filename/library metadata only, no file access'),
                    mode='scone_submission_identity')
    try:
        prior = json.loads(inputs['identity'].read_text()); assert prior['status'] == 'complete'
        for record in prior['output']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        with inputs['fields'].open() as handle:
            fields_available = {row['columnId'] for row in csv.DictReader(handle, delimiter='\t')}
        fields = ['run_accession','sample_accession','sample_alias','sample_title','sample_description',
            'experiment_accession','experiment_alias','library_name','run_alias','submitted_ftp',
            'submitted_format','submitted_md5','submitted_bytes']
        assert set(fields) <= fields_available
        url = 'https://www.ebi.ac.uk/ena/portal/api/search?'+urlencode(dict(result='read_run',
            query='study_accession="PRJNA768428"',fields=','.join(fields),format='tsv',limit=0))
        with urlopen(Request(url, headers={'User-Agent':'HaploDose-submission-identity'}), timeout=90) as response:
            raw = response.read(8388609); resolved = response.url
        assert len(raw) <= 8388608
        path = args.out/'submission_metadata.tsv'; path.write_bytes(raw)
        rows = list(csv.DictReader(io.StringIO(raw.decode()), delimiter='\t'))
        assert rows and set(fields) <= set(rows[0])
        by_run = {row['run_accession']:row for row in rows}; assert len(by_run) == len(rows)
        with inputs['earlier_runs'].open() as handle:
            earlier = {row['run_accession']:row for row in csv.DictReader(handle, delimiter='\t')}
        selected = json.loads((args.task/'identity_v1/selected.json').read_text())
        targets = []
        for gsm,srr,ordinal in selected:
            row = by_run[srr]
            assert row['sample_alias'] == gsm and row['sample_accession'] == earlier[srr]['sample_accession']
            assert row['experiment_accession'] == earlier[srr]['experiment_accession']
            targets.append(dict(row,title_number_literal=ordinal))
        (args.out/'selected_metadata.json').write_text(json.dumps(targets, indent=2)+'\n')
        (args.out/'report.json').write_text(json.dumps(dict(requested_url=url,resolved_url=resolved,
            archived_metadata=file_record(path), study_runs=len(rows), selected_runs_exact=len(targets),
            nonempty_selected_fields={key:sum(bool(row[key]) for row in targets) for key in fields},
            correspondence_rule_adopted=False, original_files_downloaded=False,
            biological_identity_certified=False, biological_null_specimens_verified=0), indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__': main()
