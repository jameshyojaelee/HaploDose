"""Inventory published DNA coordinates before assigning an interval/truth contract."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

from haplodose.provenance import start_run, finish_run, file_record


def extent_checks(start, end, length):
    return dict(zero_based_half_open=0 <= start < end <= length,
                one_based_closed=1 <= start <= end <= length)


def fixtures():
    assert extent_checks(0,10,10) == dict(zero_based_half_open=True,one_based_closed=False)
    assert extent_checks(1,10,10) == dict(zero_based_half_open=True,one_based_closed=True)
    assert extent_checks(1,11,10) == dict(zero_based_half_open=False,one_based_closed=False)
    assert extent_checks(10,10,10) == dict(zero_based_half_open=False,one_based_closed=True)
    return 4


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    root = args.task.parent
    identity = root/'scone_cell_identity_20261109/identity_v1'
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        sources=args.task/'source_records.json', headers=identity/'headers.json',
        prior_headers=identity/'run_manifest.json',
        DNA_index=root/'scone_count_chain_20261109/dna_index_v2/run_manifest.json')
    run = start_run(args.out, inputs, dict(assay_values_parsed=False),
        dict(scope='Published reference-coordinate gate only; no CNA truth or count comparison'),
        mode='scone_dna_coordinates')
    try:
        checked = fixtures()
        for key in ('prior_headers','DNA_index'):
            prior = json.loads(inputs[key].read_text()); assert prior['status'] == 'complete'
            for record in prior['output']['files'].values():
                assert file_record(record['path'])['sha256'] == record['sha256']
        sources = json.loads(inputs['sources'].read_text())['records']
        for source in sources:
            content = source['content'].encode()
            actual = hashlib.sha1(b'blob '+str(len(content)).encode()+b'\0'+content).hexdigest()
            assert actual == source['git_blob_sha'], source['path']
        header = next(record for record in json.loads(inputs['headers'].read_text())
            if '_Cellline_scWGS_' in record['requested_url'])
        assert header['status'] == 'header_available'
        url = header['requested_url']; archive = args.out/'published_DNA.tsv.gz'
        with urlopen(Request(url, headers={'User-Agent':'HaploDose-coordinate-contract'}), timeout=120) as response:
            with archive.open('xb') as handle:
                size = 0
                while block := response.read(1024*1024):
                    size += len(block); assert size <= 33554432
                    handle.write(block)
            resolved = response.url
        (args.out/'retrieval.json').write_text(json.dumps(dict(requested_url=url,resolved_url=resolved,
            archive=file_record(archive)), indent=2)+'\n')
        index = json.loads((root/'scone_count_chain_20261109/dna_index_v2/report.json').read_text())
        fai = Path(index['private_reference']['path']+'.fai')
        lengths = {row.split('\t')[0]:int(row.split('\t')[1]) for row in fai.read_text().splitlines()}
        by_chr = defaultdict(list); features = Counter(); intervals = Counter(); rows = []
        with gzip.open(archive, 'rb') as handle:
            raw_header = handle.readline(1048577)
            assert raw_header == Path(header['file']['path']).read_bytes()
            assert header['labels'][:4] == ['CHR','START','END','feature']
            assert 'One.HCT.01_dedup' in header['labels']
            for number, raw in enumerate(handle, 1):
                assert len(raw) <= 1048576 and raw.endswith(b'\n')
                fields = raw.decode('utf-8').rstrip('\r\n').split('\t',4)
                assert len(fields) == 5, 'Missing coordinate or assay fields'
                chrom, start, end, feature = fields[:4]
                assert chrom and feature
                start, end = int(start), int(end)
                rows.append((number,chrom,start,end,feature))
                by_chr[chrom].append((number,start,end)); features[feature] += 1
                intervals[chrom,start,end] += 1
        assert rows
        with (args.out/'coordinates.tsv').open('x') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['source_row','CHR_literal','START_literal','END_literal','feature_literal'])
            writer.writerows(rows)
        summary = {}; invalid = []; contiguous_blocks = Counter(); previous_chr = None
        for number,chrom,start,end,feature in rows:
            if chrom != previous_chr: contiguous_blocks[chrom] += 1
            previous_chr = chrom
        for chrom, records in by_chr.items():
            # Naming correspondence only, explicitly distinct from assembly or
            # sequence identity of the original published reference.
            candidate = chrom[3:] if chrom.startswith('chr') else chrom
            if candidate == 'M': candidate = 'MT'
            known_length = lengths.get(candidate)
            checks = Counter(); widths = Counter(); deltas = Counter()
            for position,(number,start,end) in enumerate(records):
                widths[end-start] += 1
                if position: deltas[start-records[position-1][2]] += 1
                if known_length is not None:
                    bounds = extent_checks(start,end,known_length)
                    for key,value in bounds.items(): checks[key+'_valid'] += int(value)
                    if not any(bounds.values()): invalid.append(dict(source_row=number,CHR=chrom,START=start,END=end,length=known_length))
            summary[chrom] = dict(rows=len(records), reference_contig_candidate=candidate, reference_length=known_length,
                first_record=records[0], last_record=records[-1], contiguous_source_blocks=contiguous_blocks[chrom],
                raw_END_minus_START_histogram=dict(sorted(widths.items())),
                raw_START_minus_previous_END_histogram=dict(sorted(deltas.items())),
                extent_checks=dict(checks))
        report = dict(fixtures_passed=checked, primary_source_files_verified=len(sources),
            source_rows=len(rows), header_exact=True, declared_pilot_label_present=True,
            coordinates_by_chromosome=summary,
            duplicate_features={key:n for key,n in features.items() if n>1},
            duplicate_intervals=[dict(CHR=k[0],START=k[1],END=k[2],n=n) for k,n in intervals.items() if n>1],
            intervals_invalid_under_both_extent_conventions=invalid,
            coordinate_convention_assigned=False, original_reference_sequence_identity='Not established',
            current_Ginkgo_code_proves_executed_published_version=False,
            compressed_archive_contains_other_assay_columns=True, assay_values_parsed=False,
            DNA_CNA_truth_constructed=False, biological_null_specimens_verified=0,
            independent_audit='Not run')
        (args.out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__': main()
