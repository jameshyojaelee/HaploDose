"""Collect registry metadata for the published colorectal serial sections."""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
from urllib.parse import urlencode, urlsplit
from urllib.request import Request, urlopen


ENA = 'https://www.ebi.ac.uk/ena/portal/api/'
GITHUB = 'https://api.github.com/repos/buenrostrolab/slide_dna_seq_analysis/'
ZENODO = 'https://zenodo.org/api/records/5553305'
STUDY = 'PRJNA768453'
SECTIONS = {
    'figure4_dna': 'human_colon_cancer_4_dna_200114_13',
    'figure4_rna': 'human_colon_cancer_4_rna_200102_06',
    'figure3_dna': 'human_colon_cancer_3_dna_191204_19',
    'figure3_scwgs_1': 'human_colon_cancer_3_scwgs_200227_01',
    'figure3_scwgs_2': 'human_colon_cancer_3_scwgs_210402_02',
}
FIELDS = (
    'study_accession', 'secondary_study_accession', 'sample_accession',
    'sample_alias', 'sample_title', 'experiment_accession', 'experiment_alias',
    'experiment_title', 'run_accession', 'run_alias', 'scientific_name',
    'library_strategy', 'library_source', 'library_layout', 'read_count',
    'base_count', 'fastq_bytes', 'fastq_md5',
)


def fetch(url, output, requests, *, accept=None):
    entry = {'url': url, 'file': output.name}
    requests.append(entry)
    headers = {'User-Agent': 'HaploDose-registry-metadata'}
    if accept:
        headers['Accept'] = accept
        entry['accept'] = accept
    with urlopen(Request(url, headers=headers), timeout=35) as response:
        if urlsplit(response.url).hostname != urlsplit(url).hostname:
            raise ValueError(f'Unexpected registry redirect: {response.url}')
        limit = 10 * 1024**2
        content = response.read(limit + 1)
    if len(content) > limit:
        raise ValueError('Registry metadata exceeds the response limit')
    output.write_bytes(content)
    entry.update(bytes=len(content), sha256=hashlib.sha256(content).hexdigest())
    return content


def normalize(label):
    return re.sub(r'[^a-z0-9]+', '_', label.lower()).strip('_')


def section_matches(rows):
    fields = ('sample_alias', 'sample_title', 'experiment_alias', 'experiment_title', 'run_alias')
    return {
        role: [dict(row, matched_fields=[field for field in fields
                    if normalize(row[field]) == name])
               for row in rows if any(normalize(row[field]) == name for field in fields)]
        for role, name in SECTIONS.items()
    }


def ena_metadata(out, requests):
    query = urlencode({'accession': STUDY, 'result': 'read_run',
                       'fields': ','.join(FIELDS), 'format': 'tsv'})
    content = fetch(ENA + 'filereport?' + query, out / 'ena.tsv', requests)
    rows = list(csv.DictReader(io.StringIO(content.decode()), delimiter='\t'))
    if not rows or not set(FIELDS) <= rows[0].keys():
        raise ValueError('ENA returned no runs or an incomplete schema')
    run_ids = [row['run_accession'] for row in rows]
    if any(not run_id for run_id in run_ids) or len(set(run_ids)) != len(run_ids):
        raise ValueError('Blank or duplicate run accession')
    if any(STUDY not in (row['study_accession'], row['secondary_study_accession']) for row in rows):
        raise ValueError('Run metadata contains a different BioProject')
    matches = section_matches(rows)
    section_runs = [row['run_accession'] for group in matches.values() for row in group]
    if len(set(section_runs)) != len(section_runs):
        raise ValueError('One run was assigned to multiple sections')
    return {'runs': len(rows), 'matches': matches,
            'human_titles': sorted({row['experiment_title'] for row in rows
                                    if row['scientific_name'] == 'Homo sapiens'})}


def github_metadata(out, requests):
    commit = json.loads(fetch(GITHUB + 'commits/main', out / 'github_commit.json', requests))
    tree_id = commit['commit']['tree']['sha']
    if not re.fullmatch(r'[0-9a-f]{40}', tree_id):
        raise ValueError('Invalid repository tree hash')
    tree = json.loads(fetch(GITHUB + f'git/trees/{tree_id}?recursive=1',
                            out / 'github_tree.json', requests))
    if tree['sha'] != tree_id or tree['truncated']:
        raise ValueError('Unmatched or incomplete repository listing')
    blobs = [entry for entry in tree['tree'] if entry['type'] == 'blob']
    return {'commit': commit['sha'], 'blobs': len(blobs),
            'section_paths': {role: [entry for entry in blobs if name in normalize(entry['path'])]
                              for role, name in SECTIONS.items()},
            'registration_paths': [entry for entry in blobs if
                                   re.search(r'registr|align|landmark|figure_?4|fig_?4', entry['path'], re.I)],
            'source_paths': [entry['path'] for entry in blobs if entry['path'].startswith('scripts/')]}


def zenodo_metadata(out, requests):
    record = json.loads(fetch(ZENODO, out / 'zenodo_record.json', requests))
    if str(record['id']) != '5553305':
        raise ValueError('Unmatched archive record')
    return {'doi': record['doi'], 'files': [{key: file[key] for key in ('key', 'size', 'checksum')}
                                          for file in record['files']]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    out = parser.parse_args().out
    out.mkdir()
    requests = []
    sources = {}
    for name, collect in [('ena', ena_metadata), ('github', github_metadata), ('zenodo', zenodo_metadata)]:
        try:
            sources[name] = {'status': 'available', **collect(out, requests)}
        except (OSError, ValueError, KeyError) as error:
            sources[name] = {'status': 'unresolved', 'error': repr(error)}
    report = {'sections_from_scp1278': SECTIONS, 'requests': requests, 'sources': sources,
              'measurement_files_opened': False, 'biological_specimens_evaluated': 0}
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    if all(source['status'] == 'unresolved' for source in sources.values()):
        raise SystemExit('All metadata sources unresolved; see report.json')


if __name__ == '__main__':
    main()
