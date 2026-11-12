"""Original SRA record metadata for the locked scONE pilot, without read access."""
import argparse
import csv
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

from haplodose.provenance import start_run, finish_run, file_record


def file_attributes(root):
    return [dict(node.attrib) for node in root.iter()
            if node.tag.rsplit('}', 1)[-1] in ('FILE', 'SRAFile')]


def file_attribute_fixture():
    root = ET.fromstring('<records><FILE filename="legacy.fastq"/>'
        '<SRAFile filename="native.fastq" supertype="Original"/>'
        '<ns:SRAFile xmlns:ns="urn:fixture" filename="namespaced.fastq"/></records>')
    expected = [dict(filename='legacy.fastq'), dict(filename='native.fastq', supertype='Original'),
                dict(filename='namespaced.fastq')]
    assert file_attributes(root) == expected
    return 3


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'SRA_DESIGN.md',
        selected=args.task.parent/'scone_raw_chain_20261109/chain_v1/selected_library.json',
        ENA=args.task.parent/'scone_feasibility_20261109/metadata_v1/ENA_runs.tsv')
    run = start_run(args.out, inputs, dict(locked_run='SRR16194948'),
        dict(scope='Original deposited SRA metadata only; no sequencing file access'), mode='scone_sra_identity')
    try:
        file_attribute_fixture()
        selected = json.loads(inputs['selected'].read_text())['record']
        assert selected['run_accession'] == 'SRR16194948'
        with inputs['ENA'].open() as handle:
            rows = [row for row in csv.DictReader(handle, delimiter='\t') if row['run_accession'] == selected['run_accession']]
        assert len(rows) == 1 and rows[0] == selected
        requests = []
        def fetch(endpoint, params, name):
            url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'+endpoint+'?'+urlencode(params)
            with urlopen(Request(url, headers={'User-Agent':'HaploDose-SRA-identity'}), timeout=90) as response:
                content = response.read(4194305); resolved = response.url
            assert len(content) <= 4194304
            path = args.out/name; path.write_bytes(content)
            requests.append(dict(requested_url=url, resolved_url=resolved, file=file_record(path)))
            (args.out/'requests.json').write_text(json.dumps(requests, indent=2)+'\n')
            return content
        search = json.loads(fetch('esearch.fcgi', dict(db='sra',term=selected['experiment_accession']+'[ACCN]',
            retmode='json',retmax=2), 'search.json'))['esearchresult']
        assert int(search['count']) == 1 and len(search['idlist']) == 1, 'Ambiguous or absent experiment UID'
        root = ET.fromstring(fetch('efetch.fcgi', dict(db='sra',id=search['idlist'][0],retmode='xml'), 'experiment.xml'))
        experiments = [node for node in root.iter('EXPERIMENT') if node.get('accession') == selected['experiment_accession']]
        runs = [node for node in root.iter('RUN') if node.get('accession') == selected['run_accession']]
        assert len(experiments) == len(runs) == 1
        biosamples = [node.text for node in root.iter('EXTERNAL_ID') if node.get('namespace') == 'BioSample']
        assert selected['sample_accession'] in biosamples, 'Locked BioSample identity not present'
        files = file_attributes(root)
        libraries = [node.text for node in root.iter('LIBRARY_NAME')]
        names = [dict(tag=node.tag,attributes=dict(node.attrib),text=node.text)
            for node in root.iter() if node.tag in ('TITLE','PRIMARY_ID','SUBMITTER_ID','EXTERNAL_ID')]
        (args.out/'report.json').write_text(json.dumps(dict(experiment=selected['experiment_accession'],
            run=selected['run_accession'], BioSample=selected['sample_accession'],
            experiment_run_BioSample_exact=True, library_names=libraries,
            file_attributes=files, deposited_names=names, correspondence_rule_adopted=False,
            sequencing_files_downloaded=False, biological_identity_certified=False,
            biological_null_specimens_verified=0), indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__': main()
