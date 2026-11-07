"""Trace one evidence-selected donor from channel pileups to phased allele counts.

Compute-node diagnostic, not CNV fitting or donor-genotype validation. All writes
go to a new output root. Large source matrices are read one at a time; optional
BAM inspection occurs only for mismatches with an existing ENA/channel mapping.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import io, sparse

from haplodose.alleles import read_alleles
from haplodose.io import new_output_dir, write_json
from haplodose.provenance import content_id, file_record, start_run, finish_run


SOURCE = Path(__file__).resolve().parents[2]
REQUIRED_VCF = ['CHROM', 'POS', 'REF', 'ALT', 'INFO']


def unique(values, label):
    index = pd.Index(values)
    if index.hasnans or not index.is_unique or any(not str(v).strip() for v in index):
        raise ValueError(f'{label} must be unique and nonmissing')
    return index


def vcf(path, sample=None):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        comments = []
        for line in handle:
            if line.startswith('##'):
                comments.append(line.rstrip())
            elif line.startswith('#CHROM'):
                names = line.lstrip('#').rstrip().split('\t')
                break
        else:
            raise ValueError(f'No VCF header: {path}')
        frame = pd.read_csv(handle, sep='\t', names=names, dtype={'CHROM': str, 'POS': np.int64},
                            keep_default_na=False)
    if not set(REQUIRED_VCF) <= set(frame):
        raise ValueError(f'Incomplete VCF columns: {path}')
    frame['CHROM'] = frame.CHROM.str.removeprefix('chr')
    frame['snp_id'] = (frame.CHROM + '_' + frame.POS.astype(str) + '_' + frame.REF + '_' + frame.ALT)
    unique(frame.snp_id, f'VCF SNP identities in {path}')
    for key in ['AD', 'DP', 'OTH']:
        values = frame.INFO.str.extract(r'(?:^|;)' + key + r'=(\d+)(?:;|$)', expand=False)
        frame[key + '_all'] = pd.to_numeric(values, errors='coerce')
    if sample is not None:
        if sample not in frame or 'FORMAT' not in frame or not frame.FORMAT.eq('GT').all():
            raise ValueError(f'Expected FORMAT=GT and the donor column {sample} in {path}')
        frame['GT'] = frame[sample]
    return frame, comments


def matrix(path):
    value = io.mmread(path)
    if value.data.size and (not np.isfinite(value.data).all() or np.any(value.data < 0)
                           or np.any(value.data != np.floor(value.data))
                           or np.any(value.data > np.iinfo(np.int32).max)):
        raise ValueError(f'Malformed or unsupported matrix counts in {path}')
    return value.astype(np.int64).tocsr()


def exact_sparse(a, b):
    return a.shape == b.shape and (a != b).nnz == 0


def load_barcodes(path):
    with Path(path).open() as handle:
        cells = [line.rstrip('\n') for line in handle]
    unique(cells, f'barcodes in {path}')
    return cells


def capture_numbat(data_root, out):
    prefix = data_root.parent / 'CRC/data/interim/numbat_env'
    script = prefix / 'bin/Rscript'
    if not script.is_file():
        raise ValueError('The existing Numbat R environment is unavailable; preprocessing assumptions cannot be bound')
    env = dict(os.environ, PATH=str(prefix / 'bin') + os.pathsep + os.environ.get('PATH', ''),
               R_LIBS_USER=str(prefix / 'lib/R/library'), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    env.pop('LD_LIBRARY_PATH', None)
    expression = ('a<-commandArgs(TRUE); n<-asNamespace("numbat"); sink(a[1]); '
                  'cat("numbat ",as.character(packageVersion("numbat")),"\\n"); '
                  'for (f in c("preprocess_allele","genotype","get_snps","make_vcf_chr")) '
                  '{cat("\\nFUNCTION ",f,"\\n");writeLines(deparse(get(f,envir=n)))};sink()')
    result = subprocess.run([str(script), '--vanilla', '-e', expression, str(out / 'numbat_functions.R')],
                            capture_output=True, text=True, timeout=120, env=env)
    (out / 'numbat_source_capture.log').write_text(result.stdout + result.stderr)
    result.check_returncode()
    captured = (out / 'numbat_functions.R').read_text()
    if captured.splitlines()[0].split() != ['numbat', '1.5.2']:
        raise ValueError('Trace reconstruction is specified for the inspected Numbat 1.5.2 preprocessing')
    package = prefix / 'lib/R/library/numbat'
    record = {str(p.relative_to(package)): file_record(p) for p in
              [package / 'DESCRIPTION', package / 'R/numbat.rdb', package / 'R/numbat.rdx']}
    record['Rscript_binary'] = file_record(script)
    return record


def channel_trace(channel, donor, root, out, mismatches):
    line_dir = donor / 'pileup' / channel
    pool_dir = root / 'results/kinker_numbat_20261102/pileup' / channel
    mapping = pd.read_csv(pool_dir / 'cell_map.tsv', sep='\t', dtype=str)
    for column in ['bam_barcode', 'cell_id']:
        unique(mapping[column], f'{channel} {column} map')
    all_barcodes = load_barcodes(pool_dir / 'cellSNP.samples.tsv')
    selected = mapping[mapping.line.eq(donor.name)].copy()
    j = pd.Index(all_barcodes).get_indexer(selected.bam_barcode)
    missing = selected.loc[j < 0, 'cell_id'].tolist()
    present = j >= 0
    selected = selected.loc[present].reset_index(drop=True)
    j = j[present]
    expected_cells = selected.cell_id.tolist()
    actual_cells = load_barcodes(line_dir / 'cellSNP.samples.tsv')
    if actual_cells != expected_cells:
        mismatches.append(dict(stage='barcode_order', channel=channel, detail='line barcode columns differ from mapped donor order'))
        raise ValueError(f'Unsafe count alignment: {channel} donor barcode order differs')
    # Verify the channel suffix and published builder barcode transform independently.
    channel_number, pool_number = channel.removeprefix('Pool').split('_')
    expected_suffix = f'-{channel_number}-{pool_number}'
    if (not selected.cell_id.str.endswith(expected_suffix).all()
            or not (selected.cell_id.str.split('-').str[0] + '-1').eq(selected.bam_barcode).all()):
        raise ValueError(f'{channel} barcode conversion contradicts the builder contract')
    full_vcf, _ = vcf(pool_dir / 'cellSNP.base.vcf')
    donor_vcf, _ = vcf(line_dir / 'cellSNP.base.vcf')
    subset = {}
    for name in ['AD', 'DP', 'OTH']:
        full = matrix(pool_dir / f'cellSNP.tag.{name}.mtx')
        if full.shape != (len(full_vcf), len(all_barcodes)):
            raise ValueError(f'{channel} {name} dimensions disagree with source identities')
        subset[name] = full[:, j]
        del full
    totals = {name: np.asarray(value.sum(axis=1)).ravel() for name, value in subset.items()}
    if np.any((subset['AD'] - subset['DP']).data > 0):
        raise ValueError(f'{channel} source AD exceeds DP for a selected donor cell/site')
    keep = totals['DP'] + totals['OTH'] >= 2
    expected_sites = full_vcf.snp_id[keep].tolist()
    if donor_vcf.snp_id.tolist() != expected_sites:
        mismatches.append(dict(stage='line_site_order', channel=channel,
                               detail='retained VCF differs from per-donor DP+OTH>=2 source selection'))
        raise ValueError(f'Unsafe SNP alignment: {channel} line VCF differs from source selection')
    for name in ['AD', 'DP', 'OTH']:
        if not np.array_equal(donor_vcf[name + '_all'].to_numpy(), totals[name][keep]):
            mismatches.append(dict(stage='pooled_VCF_total', channel=channel, detail=name))
    actual = {name: matrix(line_dir / f'cellSNP.tag.{name}.mtx') for name in ['AD', 'DP']}
    for name, value in actual.items():
        if not exact_sparse(value, subset[name][keep]):
            detail = dict(stage='line_count_matrix', channel=channel, detail=name)
            if value.shape == subset[name][keep].shape:
                difference = (value - subset[name][keep]).tocoo()
                if difference.nnz:
                    r, c = int(difference.row[0]), int(difference.col[0])
                    detail.update(cell=actual_cells[c], snp_id=donor_vcf.snp_id.iloc[r],
                                  observed=int(value[r, c]), expected=int(subset[name][keep][r, c]))
            mismatches.append(detail)
    pd.DataFrame(dict(column=np.arange(len(actual_cells)), cell_id=actual_cells,
                      bam_barcode=selected.bam_barcode)).to_csv(out / f'{channel}_column_order.tsv', sep='\t', index=False)
    return dict(channel=channel, source_cells=len(all_barcodes), selected_cells=len(actual_cells),
                selected_barcodes_absent_from_source=missing, source_sites=len(full_vcf),
                sites_excluded_DP_plus_OTH_below2=int((~keep).sum()), retained_sites=int(keep.sum())), donor_vcf, actual, actual_cells


def expected_final(channel, sites, counts, cells, phased):
    dp = counts['DP'].tocoo()
    if counts['AD'].shape != counts['DP'].shape or counts['DP'].shape != (len(sites), len(cells)):
        raise ValueError(f'{channel} final pileup matrix dimensions disagree')
    gt = phased.set_index('snp_id').GT.reindex(sites.snp_id)
    rate_filter = (sites.DP_all > 1) & (sites.OTH_all == 0)
    het = gt.isin(['0|1', '1|0']).to_numpy()
    positive = dp.data > 0
    # Match Matrix::summary(DP), including an explicitly stored zero if one exists.
    retained = rate_filter.to_numpy()[dp.row] & het[dp.row]
    ad = np.asarray(counts['AD'][dp.row[retained], dp.col[retained]]).ravel()
    if np.any(ad > dp.data[retained]):
        raise ValueError(f'{channel} AD exceeds DP in selected pileup')
    frame = pd.DataFrame(dict(cell=np.asarray(cells)[dp.col[retained]],
                              snp_id=sites.snp_id.to_numpy()[dp.row[retained]], AD=ad,
                              DP=dp.data[retained], GT=gt.to_numpy()[dp.row[retained]], channel=channel))
    # Numbat's distinct() removes identical rows, not conflicting measurements.
    frame = frame.drop_duplicates()
    if frame.duplicated(['cell', 'snp_id']).any():
        raise ValueError('Ambiguous duplicate cell/SNP measurements in expected final table')
    counts_out = dict(stored_DP_cell_sites=int(len(dp.data)), nonzero_DP_cell_sites=int(positive.sum()),
                      explicit_zero_DP_cell_sites=int((~positive).sum()),
                      excluded_pooled_DP_or_OTH=int((~rate_filter.to_numpy()[dp.row]).sum()),
                      excluded_no_phased_het=int((rate_filter.to_numpy()[dp.row] & ~het[dp.row]).sum()),
                      expected_final_rows=len(frame))
    if counts_out['stored_DP_cell_sites'] != sum(counts_out[k] for k in
            ['excluded_pooled_DP_or_OTH', 'excluded_no_phased_het', 'expected_final_rows']):
        raise ValueError('Cell/site inclusion and exclusions do not reconcile')
    return frame, counts_out


def compare_final(final, expected, phased, mismatches):
    required = ['cell', 'snp_id', 'CHROM', 'POS', 'REF', 'ALT', 'AD', 'DP', 'GT']
    if not set(required) <= set(final):
        raise ValueError('Final allele table lacks required measurement columns')
    if final.duplicated(['cell', 'snp_id']).any():
        raise ValueError('Duplicate final cell/SNP rows require explicit upstream reconciliation')
    for name in ['AD', 'DP', 'POS']:
        value = pd.to_numeric(final[name], errors='raise').to_numpy(float)
        if not np.isfinite(value).all() or np.any(value < 0) or np.any(value != np.floor(value)):
            raise ValueError(f'Malformed final {name}')
    key = ['cell', 'snp_id']
    joined = expected.merge(final[required], on=key, how='outer', indicator=True, suffixes=('_expected', '_final'), validate='one_to_one')
    checks = []
    for name in ['AD', 'DP', 'GT']:
        checks.append(joined[name + '_expected'].eq(joined[name + '_final']))
    correct = joined['_merge'].eq('both') & np.logical_and.reduce([x.to_numpy() for x in checks])
    for row in joined.loc[~correct].head(100).to_dict('records'):
        mismatches.append(dict(stage='final_count_or_GT', channel=row.get('channel', ''),
                               cell=row['cell'], snp_id=row['snp_id'], detail=json.dumps(row, default=str)))
    info = phased.set_index('snp_id')
    identities = final.snp_id.map(info.index.to_series()).notna()
    if not identities.all():
        mismatches.append(dict(stage='final_phase_identity', detail='final SNP absent from phased VCF'))
    present = final[identities].copy()
    reference = info.reindex(present.snp_id)
    if (not np.array_equal(present.POS.to_numpy(), reference.POS.to_numpy())
            or not np.array_equal(present.CHROM.astype(str).str.removeprefix('chr').to_numpy(), reference.CHROM.to_numpy())
            or not np.array_equal(present.REF.to_numpy(), reference.REF.to_numpy())
            or not np.array_equal(present.ALT.to_numpy(), reference.ALT.to_numpy())):
        mismatches.append(dict(stage='final_variant_identity', detail='CHROM/POS/REF/ALT differs from phased source'))
    return dict(final_rows=len(final), expected_rows=len(expected), correct_rows=int(correct.sum()),
                unexpected_final_rows=int(joined['_merge'].eq('right_only').sum()),
                missing_final_rows=int(joined['_merge'].eq('left_only').sum()),
                count_or_GT_mismatches=int((joined['_merge'].eq('both') & ~correct).sum()))


def bam_metadata(root, channels):
    manifest = SOURCE / 'benchmarks/kinker/alleles/ena_bams.tsv'
    if not manifest.is_file():
        return {}, None
    ena = pd.read_csv(manifest, sep='\t', dtype=str)
    associations = {}
    for channel in channels:
        rows = ena[ena.sample_title.eq(channel)]
        if len(rows) != 1:
            associations[channel] = dict(status='unknown_nonunique_or_missing_ENA_title')
            continue
        row = rows.iloc[0]
        bam = root / 'data/external/kinker2020_bams' / (row.run_accession + '.bam')
        sizes = dict(zip(row.submitted_ftp.split(';'), map(int, row.submitted_bytes.split(';'))))
        expected_size = next((size for name, size in sizes.items() if name.endswith('.bam')), None)
        associations[channel] = dict(status='metadata_supported_channel_mapping', run=row.run_accession,
            bam=str(bam), bam_exists=bam.is_file(), index_exists=Path(str(bam) + '.bai').is_file(),
            completion_sentinels=all(Path(str(bam) + suffix).is_file() for suffix in ['.ok', '.bai.ok']),
            expected_bytes=expected_size, observed_bytes=bam.stat().st_size if bam.is_file() else None,
            declared_md5=row.submitted_md5, full_BAM_digest='not_recomputed_by_bounded_trace',
            limit='ENA title/channel correspondence does not independently validate donor genotype or coordinate assembly')
    return associations, manifest


def bounded_forensics(mismatches, associations, root, out, limit):
    samtools = root.parent / 'CRC/data/interim/numbat_env/bin/samtools'
    records = []
    env = dict(os.environ)
    env.pop('LD_LIBRARY_PATH', None)
    for mismatch in mismatches:
        if len(records) >= limit or not samtools.is_file():
            break
        association = associations.get(mismatch.get('channel'), {})
        snp = mismatch.get('snp_id', '').split('_')
        if (len(snp) != 4 or not association.get('bam_exists') or not association.get('index_exists')
                or not association.get('completion_sentinels')
                or association.get('expected_bytes') != association.get('observed_bytes')):
            continue
        chrom, position = snp[0], int(snp[1])
        header_command = [str(samtools), 'view', '-H', association['bam']]
        try:
            header = subprocess.run(header_command, capture_output=True,
                                    text=True, timeout=30, check=True, env=env)
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as error:
            records.append(dict(command=header_command, status='bounded_header_inspection_failed', error=str(error)))
            continue
        contigs = {part[3:] for line in header.stdout.splitlines() if line.startswith('@SQ')
                   for part in line.split('\t') if part.startswith('SN:')}
        chrom = chrom if chrom in contigs else 'chr' + chrom
        if chrom not in contigs:
            continue
        region = f'{chrom}:{position}-{position}'
        mapping = pd.read_csv(root / 'results/kinker_numbat_20261102/pileup' /
                              mismatch['channel'] / 'cell_map.tsv', sep='\t', dtype=str)
        barcode = mapping.loc[mapping.cell_id.eq(mismatch.get('cell')), 'bam_barcode']
        if len(barcode) != 1:
            continue
        command = [str(samtools), 'view', '-c', '-d', 'CB:' + barcode.iloc[0], association['bam'], region]
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=30, env=env)
        except subprocess.TimeoutExpired:
            records.append(dict(command=command, status='bounded_locus_inspection_timeout'))
            continue
        records.append(dict(command=command, exit=result.returncode, stdout=result.stdout, stderr=result.stderr,
                            limit='single existing indexed-BAM locus and mapped CB tag; not UMI repileup or donor-genotype validation'))
    (out / 'bounded_bam_inspection.json').write_text(json.dumps(records, indent=2) + '\n')
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--evidence-panel', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--max-forensic-sites', type=int, default=2)
    args = parser.parse_args(argv)
    if not 0 <= args.max_forensic_sites <= 4:
        parser.error('bounded forensic site limit must be between zero and four')
    root, panel_path = args.data_root.resolve(), args.evidence_panel.resolve()
    selected = [p for p in json.loads(panel_path.read_text()) if p.get('dataset') == 'Kinker alleles']
    if len(selected) != 1 or selected[0].get('status') != 'structural_artifact_completeness_only':
        raise ValueError('Evidence panel does not identify exactly one structurally complete donor')
    selection = selected[0]
    donor = Path(selection['location']).resolve()
    if not donor.is_relative_to(root / 'results/kinker_numbat_20261102') or donor.name != selection['selected_line']:
        raise ValueError('Donor path/identity does not match the bounded diagnostic root')
    channels = sorted(p.name for p in (donor / 'pileup').iterdir() if p.is_dir())
    out_path = args.out.resolve()
    if any(out_path.is_relative_to(p) for p in [root / 'results', root / 'data', SOURCE, root.parent / 'CRC']):
        raise ValueError('Use a new diagnostic output outside source, historical results and source data')
    out = new_output_dir(out_path)
    evidence_record_path = panel_path.with_name('input_record.json')
    evidence_record = json.loads(evidence_record_path.read_text())
    if evidence_record['id'] != content_id({k: v for k, v in evidence_record.items() if k != 'id'}):
        raise ValueError('Evidence input record content ID is invalid')
    inputs = {'evidence_panel': panel_path, 'evidence_input_record': evidence_record_path, 'runner': Path(__file__),
              'builder': SOURCE / 'benchmarks/kinker/alleles/phase_pool.R',
              'pileup_builder': SOURCE / 'benchmarks/kinker/alleles/pileup_channel.sbatch'}
    for channel in channels:
        for directory in [donor / 'pileup' / channel, root / 'results/kinker_numbat_20261102/pileup' / channel]:
            for path in sorted(directory.iterdir()):
                if path.is_file() and path.name in ['cell_map.tsv', 'cellSNP.samples.tsv', 'cellSNP.base.vcf',
                                                     'cellSNP.tag.AD.mtx', 'cellSNP.tag.DP.mtx', 'cellSNP.tag.OTH.mtx']:
                    inputs[str(path)] = path
    inputs['final'] = donor / f'{donor.name}_allele_counts.tsv.gz'
    for path in sorted((donor / 'phasing').glob('*.vcf.gz')):
        inputs[str(path)] = path
    associations, ena_path = bam_metadata(root, channels)
    if ena_path is not None:
        inputs['ENA_channel_manifest'] = ena_path
    record = start_run(out, inputs, dict(max_forensic_sites=args.max_forensic_sites, selected_line=donor.name),
                       dict(measurement='existing ALT/total read counts; one-based SNP coordinates',
                            biological_unit='one selected cell line, two technical channels',
                            donor_genotype_truth='unavailable'), mode='allele_integrity_trace')
    mismatches, channel_records, expected_parts, site_parts = [], [], [], []
    try:
        current_files = {value['path']: value['sha256'] for value in record['input']['files'].values()}
        for path in selection['files']:
            expected_digest = evidence_record['files'][str(Path(path).resolve())]['sha256']
            if current_files.get(str(Path(path).resolve())) != expected_digest:
                raise ValueError(f'Allele artifact changed since evidence selection: {path}')
        numbat_source = capture_numbat(root, out)
        phased_parts, source_parts = [], []
        header_warnings = []
        for chromosome in range(1, 23):
            path = donor / 'phasing' / f'{donor.name}_chr{chromosome}.phased.vcf.gz'
            phased, header = vcf(path, donor.name)
            phased_parts.append(phased)
            unphased, _ = vcf(donor / 'phasing' / f'{donor.name}_chr{chromosome}.vcf.gz', donor.name)
            source_parts.append(unphased)
            if not set(phased.snp_id) <= set(unphased.snp_id):
                mismatches.append(dict(stage='phase_source', detail=f'chr{chromosome}: phased SNP absent from genotype input'))
            if chromosome == 1:
                header_warnings = [line for line in header if line.startswith('##contig=<ID=chr1,')]
        phased = pd.concat(phased_parts, ignore_index=True)
        genotype_source = pd.concat(source_parts, ignore_index=True)
        unique(phased.snp_id, 'phased SNP identities')
        for channel in channels:
            stats, sites, counts, cells = channel_trace(channel, donor, root, out, mismatches)
            expected, reconciliation = expected_final(channel, sites, counts, cells, phased)
            channel_records.append(dict(**stats, **reconciliation))
            expected_parts.append(expected)
            site_parts.append(sites)
        totals = pd.concat(site_parts).groupby('snp_id')[['AD_all', 'DP_all', 'OTH_all']].sum()
        ratio = totals.AD_all / totals.DP_all
        eligible = ratio.between(.1, .9) | ((ratio == 1) & (totals.DP_all >= 10))
        expected_genotypes = pd.Series(np.where(ratio.between(.1, .9), '0/1', '1/1'), index=totals.index)[eligible]
        observed_genotypes = genotype_source.set_index('snp_id').GT
        if not expected_genotypes.sort_index().equals(observed_genotypes.sort_index()):
            mismatches.append(dict(stage='genotyping_source_selection', detail='pooled per-donor genotype inclusion or GT differs from inspected Numbat rule'))
        expected = pd.concat(expected_parts, ignore_index=True)
        final = pd.read_csv(inputs['final'], sep='\t', dtype={'cell': str, 'CHROM': str, 'GT': str})
        reconciliation = compare_final(final, expected, phased, mismatches)
        final[['cell', 'snp_id']].to_csv(out / 'final_row_order.tsv.gz', sep='\t', index=False)
        # Independent orientation identity: HaploDose's B channel is ALT for 0|1, REF for 1|0.
        cells = unique(final.cell.drop_duplicates(), 'final cell identities')
        allele = read_alleles(inputs['final'], cells)
        site_keys = 'chr' + final.CHROM.str.removeprefix('chr') + ':' + final.POS.astype(str)
        oriented = final.assign(site=site_keys, B=np.where(final.GT.eq('0|1'), final.AD, final.DP - final.AD))
        oriented = oriented[oriented.DP > 0]
        pooled = oriented.groupby(['cell', 'site'])[['B', 'DP']].sum()
        source_sites = pd.Index(allele.snps.chromosome + ':' + allele.snps.pos.astype(str))
        row = pd.Index(cells).get_indexer(pooled.index.get_level_values('cell'))
        column = source_sites.get_indexer(pooled.index.get_level_values('site'))
        if np.any(row < 0) or np.any(column < 0):
            raise ValueError('Oriented counts lost a final cell/site identity')
        if (not np.array_equal(np.asarray(allele.A[row, column]).ravel(), pooled.B.to_numpy())
                or not np.array_equal(np.asarray(allele.D[row, column]).ravel(), pooled.DP.to_numpy())):
            mismatches.append(dict(stage='HaploDose_GT_orientation', detail='B or total channel disagrees with direct per-GT orientation'))
        pd.DataFrame(mismatches).to_csv(out / 'mismatches.tsv', sep='\t', index=False)
        pd.DataFrame(channel_records).to_csv(out / 'channel_reconciliation.tsv', sep='\t', index=False)
        forensic = bounded_forensics(mismatches, associations, root, out, args.max_forensic_sites) if mismatches else []
        summary = dict(selected_line=donor.name, channels=channel_records, final_reconciliation=reconciliation,
            genotype_source_sites=len(genotype_source), phased_sites=len(phased),
            phased_source_sites_dropped=len(set(genotype_source.snp_id) - set(phased.snp_id)),
            mismatches=len(mismatches), passed=not mismatches, source_bam_associations=associations,
            forensic_checks=len(forensic), installed_numbat_source=numbat_source,
            evidence_input_id=evidence_record['id'], haplodose_reader_exclusions=allele.excluded,
            orientation_checked_positive_DP_rows=len(oriented),
            phase_header_contig_chr1=header_warnings,
            coordinate_limit='VCF contig headers were inherited by Numbat; declared hg19 recipe is not independent assembly validation',
            limits=['artifact consistency does not validate donor genotype, phase accuracy, ASE bias or DNA clone identity',
                    'gene and cM annotation accuracy is outside this count/orientation trace',
                    'one selected cell line; channels/cells/SNPs are not independent biological specimens'])
        write_json(summary, out / 'summary.json')
        if mismatches:
            raise ValueError(f'{len(mismatches)} allele integrity mismatches; inspect summary and mismatch table')
        finish_run(out, record, backend='cpu')
    except Exception as error:
        pd.DataFrame(mismatches).to_csv(out / 'mismatches.tsv', sep='\t', index=False)
        finish_run(out, record, error=error)
        raise
    print(json.dumps(dict(selected_line=donor.name, passed=True, final_rows=len(final))))


if __name__ == '__main__':
    main()
