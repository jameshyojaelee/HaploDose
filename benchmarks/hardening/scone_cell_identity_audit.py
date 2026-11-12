"""Independent literal scONE metadata audit and explicitly provisional label crosswalk."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path

from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    root = args.task.parent
    inputs = dict(script=Path(__file__), design=args.task/'AUDIT_DESIGN.md',
        identity=args.task/'identity_v1/run_manifest.json',
        original_pilot=root/'scone_raw_chain_20261109/chain_v1/run_manifest.json',
        ENA=root/'scone_feasibility_20261109/metadata_v1/ENA_runs.tsv',
        author_source=root/'scone_feasibility_20261109/source_records.json')
    run = start_run(args.out, inputs, {}, dict(scope='Administrative provenance, not biological correspondence'),
                    mode='scone_cell_identity_audit')
    try:
        identity = json.loads(inputs['identity'].read_text()); assert identity['status'] == 'complete'
        for record in identity['output']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        pilot = json.loads(inputs['original_pilot'].read_text()); assert pilot['status'] == 'complete'
        pilot_path = root/'scone_raw_chain_20261109/chain_v1/sample.soft.txt'
        earlier = [r for r in pilot['output']['files'].values() if r['path'] == str(pilot_path)]
        assert len(earlier) == 1 and file_record(pilot_path)['sha256'] == earlier[0]['sha256']
        with inputs['ENA'].open() as handle:
            ena = list(csv.DictReader(handle, delimiter='\t'))
        selected = {}
        for row in ena:
            name = row['experiment_title'].split(': ',2)[-1]
            prefix, suffix = 'HCT116 scONEseq ', ' Homo sapiens OTHER'
            if name.startswith(prefix) and name.endswith(suffix):
                number = name[len(prefix):-len(suffix)]
                assert number.isdecimal() and row['sample_alias'] not in selected
                selected[row['sample_alias']] = (row['run_accession'],number)
        saved = json.loads((args.task/'identity_v1/selected.json').read_text())
        assert [[gsm,*value] for gsm,value in sorted(selected.items())] == saved
        summaries = {row['GSM']:row for row in json.loads((args.task/'identity_v1/samples.json').read_text())}
        requests = json.loads((args.task/'identity_v1/requests.json').read_text())
        assert len(requests) == len(selected) and len({row['GSM'] for row in requests}) == len(requests)
        source_records = {}; offsets = Counter(); description_ids = Counter()
        for request in requests:
            assert request['GSM'] in selected
            if request['status'] != 'metadata_available': continue
            gsm = request['GSM']; path = Path(request['file']['path'])
            assert file_record(path)['sha256'] == request['file']['sha256']
            fields = {}
            for line in path.read_text().splitlines():
                left, separator, right = line.partition(' = ')
                if separator: fields.setdefault(left, []).append(right)
            assert fields['^SAMPLE'] == [gsm]
            run_id, ordinal = selected[gsm]
            assert fields['!Sample_title'] == ['HCT116 scONEseq '+ordinal]
            descriptions = fields.get('!Sample_description', [])
            numbers = [value[len('HCT116 single cell '):] for value in descriptions
                if value.startswith('HCT116 single cell ') and value[len('HCT116 single cell '):].isdecimal()]
            saved_row = summaries[gsm]
            assert saved_row['SRR'] == run_id and saved_row['title_number_literal'] == ordinal
            assert saved_row['description_numbers_literal'] == numbers
            assert saved_row['descriptions'] == descriptions
            assert saved_row['source_names'] == fields.get('!Sample_source_name_ch1', [])
            assert saved_row['relations'] == fields.get('!Sample_relation', [])
            if len(numbers) == 1:
                offsets[int(numbers[0])-int(ordinal)] += 1; description_ids[numbers[0]] += 1
            source_records[gsm] = (run_id,ordinal,numbers)
        assert set(source_records) == set(summaries)
        headers = json.loads((args.task/'identity_v1/headers.json').read_text())
        label_sets = {}
        for record in headers:
            if record['status'] != 'header_available': continue
            path = Path(record['file']['path']); assert file_record(path)['sha256'] == record['file']['sha256']
            with path.open(encoding='utf-8-sig', newline='') as handle:
                labels = list(csv.reader(handle, delimiter=record['delimiter']))
            assert len(labels) == 1 and labels[0] == record['labels']
            assert {key:n for key,n in Counter(labels[0]).items() if n>1} == record['duplicate_labels']
            kind = 'RNA' if '_RNA_' in record['requested_url'] else 'DNA'
            assert kind not in label_sets
            label_sets[kind] = set(labels[0])
        assert set(label_sets) == {'RNA','DNA'}
        # These are explicit metadata-only hypotheses, not adopted identities.
        candidate_rows = []; supports = Counter()
        for gsm,(run_id,title,descriptions) in sorted(source_records.items()):
            for origin,numbers in [('title',[title]),('description',descriptions)]:
                if len(numbers) != 1: continue
                label = 'One.HCT.'+str(int(numbers[0])).zfill(2)
                rna_present = label in label_sets['RNA']; dna_present = label+'_dedup' in label_sets['DNA']
                supports[origin+'_RNA_labels_present'] += int(rna_present)
                supports[origin+'_DNA_labels_present'] += int(dna_present)
                candidate_rows.append([gsm,run_id,origin,numbers[0],label,label+'_dedup',int(rna_present),int(dna_present),
                    'provisional; biological identity unverified'])
        with (args.out/'candidate_crosswalk.tsv').open('x') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['GSM','SRR','ordinal_source','literal_ordinal','RNA_label_hypothesis','DNA_label_hypothesis',
                'RNA_header_present','DNA_header_present','status'])
            writer.writerows(candidate_rows)
        earlier_report = json.loads((args.task/'identity_v1/report.json').read_text())
        assert earlier_report['metadata_records_available'] == len(source_records)
        assert {int(key):n for key,n in earlier_report['description_title_offsets_descriptive'].items()} == dict(offsets)
        assert earlier_report['unique_description_numbers'] == len(description_ids)
        assert earlier_report['all_selected_metadata_available'] == (len(source_records) == len(selected))
        assert earlier_report['headers_available'] == len(label_sets)
        (args.out/'report.json').write_text(json.dumps(dict(selected_records_exact=len(selected),
            sample_records_exact=len(source_records), matrix_headers_exact=len(label_sets),
            descriptive_offsets=dict(offsets), duplicated_description_ordinals={k:n for k,n in description_ids.items() if n>1},
            candidate_header_support=dict(supports),
            correspondence_rule_adopted=False, expression_DNA_values_parsed=False,
            biological_identity_certified=False, biological_null_specimens_verified=0,
            decision='Metadata favors title convention; exact author crosswalk or independent molecular correspondence still needed'), indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__': main()
