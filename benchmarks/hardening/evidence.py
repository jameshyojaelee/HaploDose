"""Inventory the fixed development panel without fitting or rescoring historical results.

Run on an I/O compute node: hashes include the actual large measurement sources.
The four records describe this collection, not a completed inference experiment.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import sys

import h5py

from haplodose import provenance


def table(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def first_line(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        return handle.readline().rstrip('\n').split('\t')


def signed(record):
    return dict(record, id=provenance.content_id(record))


def select_alleles(root):
    """Lexical selection by line, then source path; existence is not semantic QC."""
    candidates = []
    base = root / 'results/kinker_numbat_20261102'
    for pool in sorted(base.glob('pool*')):
        for line in sorted((pool / 'alleles').glob('*')):
            if not line.is_dir():
                continue
            final = line / f'{line.name}_allele_counts.tsv.gz'
            phased = [line / 'phasing' / f'{line.name}_chr{c}.phased.vcf.gz' for c in range(1, 23)]
            pileups = sorted((line / 'pileup').glob('*'))
            required = [p / f for p in pileups for f in
                        ['cellSNP.base.vcf', 'cellSNP.samples.tsv', 'cellSNP.tag.AD.mtx', 'cellSNP.tag.DP.mtx']]
            files = [final, *phased, *required]
            if pileups and all(p.is_file() and p.stat().st_size > 0 for p in files):
                candidates.append((line.name, str(line), files))
    if not candidates:
        return dict(status='no_structurally_complete_autosomal_artifacts', raw_bam_association='unknown'), []
    name, location, files = min(candidates, key=lambda x: (x[0], x[1]))
    return dict(status='structural_artifact_completeness_only', selected_line=name, location=location,
                selection='lexically first line, tie broken by full artifact directory',
                phase_chromosomes=list(range(1, 23)), allele_header=first_line(files[0]),
                raw_bam_association='unknown; existing BAM files are not evidence of this donor/channel association',
                biological_unit='cell line; chromosomes and cells are not independent donors'), files


class Inventory:
    def __init__(self):
        self.files = {}

    def add(self, path, purpose):
        path = Path(path).resolve()
        key = str(path)
        if key not in self.files:
            self.files[key] = (dict(provenance.file_record(path), purposes=[]) if path.is_file()
                               else dict(path=key, status='missing', purposes=[]))
        if purpose not in self.files[key]['purposes']:
            self.files[key]['purposes'].append(purpose)
        return key


def matrix_header(path):
    with h5py.File(path, 'r') as handle:
        group = handle['matrix'] if 'matrix' in handle else handle
        cells = group['barcodes' if 'matrix' in handle else 'cell_ids'].asstr()[:].tolist()
        genes = group['features/name' if 'matrix' in handle else 'genes'].asstr()[:].tolist()
        counts = group if 'matrix' in handle else group['X']
        return dict(cell_ids_in_matrix_order=cells, feature_ids_in_matrix_order=genes,
                    sparse_data_dtype=str(counts['data'].dtype), sparse_nnz=len(counts['data']),
                    feature_filter='not executed; these are input feature identities, not fitted features')


def bulk_panel(root, directory, ids, inventory, *, counts_kind, include_all_normals=True):
    files = {name: inventory.add(directory / name, 'prepared bulk input') for name in
             ['counts.tsv', 'samples.tsv', 'positions.tsv']}
    if not (directory / 'samples.tsv').is_file() or not (directory / 'counts.tsv').is_file():
        return dict(status='missing_prepared_input', files=files, requested_ids=ids, counts_kind=counts_kind)
    rows, header = table(directory / 'samples.tsv'), first_line(directory / 'counts.tsv')
    ordered = [row for row in rows if row['sample'] in ids or (include_all_normals and row['role'] == 'normal')]
    return dict(status='inventoried', files=files, requested_ids=ids, counts_kind=counts_kind,
                count_sample_order=header[1:], selected_metadata_in_source_order=ordered,
                missing_requested_ids=sorted(set(ids) - {r['sample'] for r in rows}),
                reference_membership=[r['sample'] for r in ordered if r['role'] == 'normal'],
                gene_namespace='symbol', genome_build='hg38', coordinate_system='one-based-closed',
                measurement_validation='builder/source declaration; no full count-model validation executed')


def collect(root, source, task, out):
    inventory = Inventory()
    checksum_path = task / 'frozen_artifact_checksums.json'
    expected = json.loads(checksum_path.read_text())

    def check_frozen():
        checked = {}
        for relative, digest in expected.items():
            actual = provenance.file_record(root / relative)
            if actual['sha256'] != digest:
                raise RuntimeError(f'Historical artifact changed: {relative}')
            checked[relative] = actual
        return checked

    frozen_before = check_frozen()
    claims_path = source / 'benchmarks/hardening/claims.json'
    claims = json.loads(claims_path.read_text())
    inventory.add(claims_path, 'frozen panel and interpretation contract')
    panels = []
    kinker = root / 'results/kinker_input_v2_20261101_1551'
    kspec = next(x for x in claims['diagnostic_panel'] if x['dataset'] == 'Kinker')
    entry = dict(dataset='Kinker', specification=kspec, counts_kind='reconstructed_integer_from_published_CPM',
                 raw_count_origin='unverified; integer reconstruction does not establish original molecule counts',
                 source_files=[inventory.add(kinker / p, 'Kinker prepared measurement/reference/truth') for p in
                               ['counts.h5', 'positions.tsv', 'depmap_cn.tsv.gz', 'eligible_lines.tsv']], rotations=[])
    if (kinker / 'counts.h5').is_file():
        entry['matrix'] = matrix_header(kinker / 'counts.h5')
    for rotation in (0, 1):
        path = kinker / f'cells_r{rotation}.tsv'
        inventory.add(path, 'Kinker matrix-order metadata')
        if not path.is_file():
            continue
        rows = table(path)
        query = [(i, r) for i, r in enumerate(rows) if r['role'] == 'query' and
                 any(r['specimen'] == name or r['specimen'].startswith(name + '_') for name in kspec['ids'])]
        datasets = {r['dataset'] for _, r in query}
        query_indices = {j for j, _ in query}
        selected = [(i, r) for i, r in enumerate(rows) if i in query_indices or
                    (r['role'] == 'reference' and r['dataset'] in datasets)]
        entry['rotations'].append(dict(rotation=rotation, selected_rows_in_matrix_order=[i for i, _ in selected],
                                       selected_metadata_in_matrix_order=[r for _, r in selected],
                                       resolved_query_ids=sorted({r['specimen'] for _, r in query})))
    panels.append(entry)

    uspec = next(x for x in claims['diagnostic_panel'] if x['dataset'] == 'U-CAN')
    panels.append(dict(dataset='U-CAN', specification=uspec, **bulk_panel(
        root, root / 'results/ucan_bulk_input_20261101_1352', uspec['ids'], inventory,
        counts_kind='raw_featureCounts_exon_counts_builder_declared')))
    tspec = next(x for x in claims['diagnostic_panel'] if x['dataset'] == 'TCGA')
    tcga = root / 'results/tcga_bulk_20261101_1405/input'
    panels.append(dict(dataset='TCGA', specification=tspec, eligible_for_count_likelihood=False, **bulk_panel(
        root, tcga / 'LUAD', tspec['ids'], inventory, counts_kind='normalized_inverted_batch_adjusted_RSEM')))
    raw, xena = tcga / 'coadread_raw', tcga / 'coadread_xena'
    matched = []
    if (raw / 'samples.tsv').is_file() and (xena / 'samples.tsv').is_file():
        matched = sorted({r['sample'] for r in table(raw / 'samples.tsv')} &
                         {r['sample'] for r in table(xena / 'samples.tsv')})
    panels.append(dict(dataset='TCGA COAD/READ', matched_sample_ids=matched,
                       partition_rule='not allocated here; fixed salted-ID partitions in frozen profile.py',
                       raw=bulk_panel(root, raw, matched, inventory, counts_kind='raw_STAR_unstranded_builder_declared',
                                      include_all_normals=False),
                       normalized=bulk_panel(root, xena, matched, inventory, counts_kind='normalized_RSEM_diagnostic_only',
                                             include_all_normals=False)))

    aspec = next(x for x in claims['diagnostic_panel'] if x['dataset'] == 'Andor')
    for line in aspec['ids']:
        directory = root / 'results/andor_20261101_1224/rna_input_v2' / line
        panel = dict(dataset='Andor', line=line, specification=aspec, counts_kind='raw_10x_UMI_builder_declared',
                     files=[inventory.add(directory / f, 'Andor prepared RNA') for f in
                            ['counts.h5', 'cells.tsv', 'positions.tsv']])
        if (directory / 'counts.h5').is_file():
            panel['matrix'] = matrix_header(directory / 'counts.h5')
        if (directory / 'cells.tsv').is_file():
            panel['selected_metadata_in_matrix_order'] = table(directory / 'cells.tsv')
        depth = 'scdna_fulldepth' if line == 'SNU-601' else 'scdna_lowdepth'
        truth = root / 'results/andor_20261101_1224' / depth / line
        panel['historical_truth'] = [inventory.add(truth / f, 'historical DNA truth, superseded diagnostic') for f in
                                     ['qc.json', 'cells.tsv', 'clone_bins.tsv']]
        panel['truth_limits'] = 'Different RNA/DNA cells; inferred DNA groups and depth-limited resolution; no cell-level RNA clone labels'
        panels.append(panel)
    allele, allele_files = select_alleles(root)
    allele['files'] = [inventory.add(path, 'selected donor allele/phasing/pileup artifact') for path in allele_files]
    if allele_files:
        allele['pileup_cell_order_by_channel'] = {}
        for path in allele_files:
            if path.name == 'cellSNP.samples.tsv':
                with path.open() as handle:
                    allele['pileup_cell_order_by_channel'][str(path.parent)] = [line.rstrip('\n') for line in handle]
    panels.append(dict(dataset='Kinker alleles', **allele))

    # Explicit transitive measurement/truth sources; no recursive scan of data or result trees.
    crc = root.parent / 'CRC'
    sources = [root / 'data/external/kinker2020_gse157220/GSE157220_CPM_data.txt.gz',
               root / 'data/external/kinker2020_gse157220/GSE157220_Pool_composition.xlsx',
               crc / 'data/external/depmap_24q4/Model.csv', crc / 'data/external/depmap_24q4/OmicsCNGene.csv',
               crc / 'data/external/depmap_24q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv',
               crc / 'data/interim/full_transcriptome_counts/exon_counts.tsv.gz',
               crc / 'data/external/ensembl112/Homo_sapiens.GRCh38.112.gtf.gz',
               crc / 'results/cnv_caller_benchmark/SMC_full/input/positions.tsv',
               crc / 'results/cnv_caller_benchmark/YCC230/input/positions.tsv',
               crc / 'data/external/tcga_pancan_xena/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz',
               crc / 'data/external/tcga_pancan_xena/TCGA_phenotype_denseDataOnlyDownload.tsv.gz',
               crc / 'data/external/tcga_pancan_cn/TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
               crc / 'data/external/tcga_pancan_cn/TCGA_mastercalls.abs_segtabs.fixed.txt',
               root / 'data/external/andor2020_gse142750/GSE142750_RAW.tar',
               root / 'data/external/andor2020_scdna/bin_annotation_500kb.tsv']
    for line in aspec['ids']:
        depth = 'bins_original_bam' if line == 'SNU-601' else 'bins'
        sources.extend(root / 'data/external/andor2020_scdna' / depth / f'{line}{suffix}' for suffix in ['.tsv.gz', '.rg.txt'])
    for sample in uspec['ids']:
        patient = sample.split('.')[2]
        sources.extend(sorted((crc / 'data/external/eva_cnv/vcf').glob(f'*CRC-SW-{patient}-T*.vcf.gz')))
    mapping = crc / 'data/external/tcga_coadread/manifests/file_case_map.tsv'
    sources.append(mapping)
    if mapping.is_file():
        for row in table(mapping):
            if row.get('data_type') == 'Gene Expression Quantification' and row.get('sample_id', '')[:15] in matched:
                sources.append(crc / 'data/external/tcga_coadread/raw/gene_expression_quantification' /
                               row['file_id'] / row['file_name'])
    for path in sources:
        inventory.add(path, 'transitive measurement/annotation/truth source')

    code = {}
    for name in ['python', 'benchmarks', 'scripts']:
        for path in sorted((source / name).rglob('*')):
            if path.is_file() and path.suffix in {'.py', '.R', '.json', '.tsv', '.sbatch'} and '__pycache__' not in path.parts:
                code[str(path.relative_to(source))] = provenance.file_record(path)
    for name in ['source_manifest.json', 'SNAPSHOT_COMPLETE', 'pyproject.toml', 'crates/core/Cargo.lock']:
        if (source / name).is_file():
            code[name] = provenance.file_record(source / name)
    input_record = signed(dict(kind='diagnostic_input_bundle', files=inventory.files, panels=panels,
                               panel_contract=claims, inference_status='not_run_by_collector'))
    execution = signed(dict(kind='evidence_collection_execution', input_id=input_record['id'], frozen_source=code,
                            actual_collector=provenance.file_record(__file__),
                            actual_provenance_helper=provenance.file_record(provenance.__file__),
                            executable=sys.executable, command=sys.argv,
                            slurm_job_id=os.environ.get('SLURM_JOB_ID'),
                            semantics='Collection/hash/header inspection only; no inference or scoring'))
    (out / 'panel.json').write_text(json.dumps(panels, indent=2) + '\n')
    output = signed(dict(kind='evidence_collection_output', execution_id=execution['id'],
                         files={'panel.json': provenance.file_record(out / 'panel.json')}, inference_status='not_run'))
    evaluation = signed(dict(kind='truth_and_scorer_contract', input_id=input_record['id'], output_id=output['id'],
                             scorer_source_id=provenance.content_id(code), status='not_scored',
                             limits=[x for x in claims['diagnostic_panel']],
                             uncertainty='No new accuracy, calibration or biological-validation claim follows from artifact collection',
                             frozen_artifacts_before=frozen_before, frozen_artifacts_after=check_frozen()))
    for name, record in [('input', input_record), ('execution', execution), ('output', output), ('evaluation', evaluation)]:
        (out / f'{name}_record.json').write_text(json.dumps(record, indent=2) + '\n')
    return dict(input_id=input_record['id'], execution_id=execution['id'], output_id=output['id'],
                evaluation_id=evaluation['id'], files=len(inventory.files), panels=len(panels))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['root', 'source', 'task-root', 'out']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args(argv)
    root, source, task, out = (p.resolve() for p in [args.root, args.source, args.task_root, args.out])
    if any(out.is_relative_to(p) for p in [root / 'results', root / 'data', source, task / 'baseline']):
        raise ValueError('Evidence output must be fresh and outside historical data/results and frozen source')
    if not (source / 'SNAPSHOT_COMPLETE').is_file():
        raise ValueError('Evidence requires a completed candidate source snapshot')
    manifest = json.loads((source / 'source_manifest.json').read_text())
    if (source / 'SNAPSHOT_COMPLETE').read_text().strip() != provenance.content_id(manifest):
        # freeze_validation.py uses the standard JSON encoding rather than content_id's compact encoding.
        import hashlib
        expected = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()
        if (source / 'SNAPSHOT_COMPLETE').read_text().strip() != expected:
            raise ValueError('Source snapshot completion marker does not match its manifest')
    for relative, digest in manifest.items():
        path = (source / relative).resolve()
        if not path.is_relative_to(source) or provenance.file_record(path)['sha256'] != digest:
            raise ValueError(f'Frozen source manifest mismatch: {relative}')
    out.mkdir(parents=True, exist_ok=False)
    try:
        summary = collect(root, source, task, out)
        (out / 'summary.json').write_text(json.dumps(dict(status='complete', **summary), indent=2) + '\n')
    except Exception as error:
        (out / 'failure.json').write_text(json.dumps(dict(status='failed', error_type=type(error).__name__,
                                                         error=str(error)), indent=2) + '\n')
        raise
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
