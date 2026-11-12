"""Independent full XML file-element audit for the locked scONE pilot."""
import argparse
import importlib.util
import json
from pathlib import Path
import xml.dom.minidom as dom

from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    helper = Path(__file__).with_name('scone_sra_identity.py')
    inputs = dict(script=Path(__file__), corrected_helper=helper, design=args.task/'SRA_AUDIT_DESIGN.md',
        original=args.task/'sra_v1/run_manifest.json', xml=args.task/'sra_v1/experiment.xml',
        headers=args.task/'identity_v1/headers.json')
    run = start_run(args.out, inputs, {}, dict(scope='Primary metadata file inventory and pilot naming link'),
                    mode='scone_sra_identity_audit')
    try:
        original = json.loads(inputs['original'].read_text()); assert original['status'] == 'complete'
        for record in original['output']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        spec = importlib.util.spec_from_file_location('corrected_sra_collector', helper)
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        assert module.file_attribute_fixture() == 3
        raw = inputs['xml'].read_bytes(); document = dom.parseString(raw)
        nodes = list(document.getElementsByTagName('*'))
        files = [dict(node.attributes.items()) for node in nodes if node.localName in ('FILE','SRAFile')]
        parsed = module.ET.fromstring(raw)
        assert module.file_attributes(parsed) == files
        experiment = [node for node in nodes if node.localName == 'EXPERIMENT' and node.getAttribute('accession') == 'SRX12479606']
        runs = [node for node in nodes if node.localName == 'RUN' and node.getAttribute('accession') == 'SRR16194948']
        assert len(experiment) == len(runs) == 1
        biosamples = [node.firstChild.data for node in nodes if node.localName == 'EXTERNAL_ID'
            and node.getAttribute('namespace') == 'BioSample' and node.firstChild]
        assert 'SAMN22044808' in biosamples
        original_fastq = [item for item in files if item.get('supertype') == 'Original' and item.get('semantic_name') == 'fastq']
        expected_names = ['One-HCT-01_read1.fastq.gz','One-HCT-01_read2.fastq.gz']
        assert sorted(item['filename'] for item in original_fastq) == expected_names
        earlier = json.loads((args.task/'sra_v1/report.json').read_text())
        old_FILE = [dict(node.attributes.items()) for node in nodes if node.localName == 'FILE']
        assert old_FILE == earlier['file_attributes']
        # This checks an explicitly declared naming hypothesis, not biological
        # identity or an automatic general rule for other deposited samples.
        headers = json.loads(inputs['headers'].read_text())
        candidates = ['One.HCT.01','One.HCT.97']; support = {}
        for record in headers:
            assert record['status'] == 'header_available'
            assert file_record(record['file']['path'])['sha256'] == record['file']['sha256']
            kind = 'RNA' if '_RNA_' in record['requested_url'] else 'DNA'
            suffix = '' if kind == 'RNA' else '_dedup'
            support[kind] = {label:label+suffix in record['labels'] for label in candidates}
        result = dict(XML_source=file_record(inputs['xml']), schema_fixture_cases=3,
            full_file_inventory_exact=files, original_fastq_metadata=original_fastq,
            earlier_FILE_only_inventory_exact_but_incomplete=True,
            experiment_run_BioSample_exact=True, candidate_header_support=support,
            raw_submitted_stem='One-HCT-01', title_based_matrix_candidate='One.HCT.01',
            naming_hypothesis='Strip explicit _read1/_read2 FASTQ suffix; replace hyphens with periods',
            pilot_administrative_naming_link_supported=True,
            biological_identity_certified=False, cohort_wide_crosswalk_adopted=False,
            sequencing_files_downloaded=False, biological_null_specimens_verified=0)
        (args.out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__': main()
