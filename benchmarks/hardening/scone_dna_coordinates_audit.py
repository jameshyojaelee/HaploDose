"""Independent literal interval/FAI audit of the published scONE DNA table."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path

from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'AUDIT_DESIGN.md',
        original=args.task/'coordinates_v1/run_manifest.json',
        sources=args.task/'interval_source_records.json',
        index=args.task.parent/'scone_count_chain_20261109/dna_index_v2/run_manifest.json')
    run = start_run(args.out, inputs, {}, dict(scope='Reference coordinates/support only; no genotype truth'),
                    mode='scone_dna_coordinates_audit')
    try:
        unique_fixture = Counter({('1',0,10):1, ('1',10,20):1})
        duplicate_fixture = Counter({('1',0,10):2})
        assert all(n == 1 for n in unique_fixture.values())
        assert not all(n == 1 for n in duplicate_fixture.values())
        for key in ('original','index'):
            original = json.loads(inputs[key].read_text()); assert original['status'] == 'complete'
            for record in original['output']['files'].values():
                assert file_record(record['path'])['sha256'] == record['sha256']
        sources = json.loads(inputs['sources'].read_text())['records']
        for source in sources:
            raw = source['content'].encode()
            assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest() == source['git_blob_sha']
        report = json.loads((args.task/'coordinates_v1/report.json').read_text())
        index = json.loads((args.task.parent/'scone_count_chain_20261109/dna_index_v2/report.json').read_text())
        fai = Path(index['private_reference']['path']+'.fai')
        lengths = {line.split('\t')[0]:int(line.split('\t')[1]) for line in fai.read_text().splitlines()}
        with (args.task/'coordinates_v1/coordinates.tsv').open(newline='') as handle:
            exported = list(csv.DictReader(handle, delimiter='\t'))
        observed = []; per_chrom = defaultdict(list); feature_counts = Counter(); interval_counts = Counter()
        with gzip.open(args.task/'coordinates_v1/published_DNA.tsv.gz', 'rt', newline='') as handle:
            reader = csv.reader(handle, delimiter='\t'); header = next(reader)
            assert header[:4] == ['CHR','START','END','feature'] and 'One.HCT.01_dedup' in header
            for number,row in enumerate(reader,1):
                # Other fields are unparsed strings, never numerical assay data.
                chrom = row[0]; first = int(row[1]); last = int(row[2]); feature = row[3]
                observed.append(dict(source_row=str(number),CHR_literal=chrom,START_literal=str(first),
                    END_literal=str(last),feature_literal=feature))
                per_chrom[chrom].append((number,first,last)); feature_counts[feature] += 1
                interval_counts[chrom,first,last] += 1
        assert observed == exported and len(observed) == report['source_rows']
        assert {key:n for key,n in feature_counts.items() if n>1} == report['duplicate_features']
        assert all(n == 1 for n in interval_counts.values()) and not report['duplicate_intervals']
        zero_based = []; summary = {}; omitted = []; width_counts = Counter()
        for chrom,rows in per_chrom.items():
            name = chrom.removeprefix('chr')
            if name == 'M': name = 'MT'
            assert name in lengths
            expected = report['coordinates_by_chromosome'][chrom]
            assert expected['rows'] == len(rows)
            assert expected['first_record'] == list(rows[0]) and expected['last_record'] == list(rows[-1])
            deltas = Counter(); widths = Counter(); valid0 = 0; valid1 = 0
            prior = None
            for number,first,last in rows:
                valid0 += int(0 <= first and first < last and last <= lengths[name])
                valid1 += int(1 <= first and first <= last and last <= lengths[name])
                widths[last-first] += 1; width_counts[last-first+1] += 1
                if prior is not None: deltas[first-prior] += 1
                prior = last
                assert first >= 1 and first <= last <= lengths[name]
                zero_based.append((number,name,first-1,last,observed[number-1]['feature_literal']))
            assert {int(k):n for k,n in expected['raw_END_minus_START_histogram'].items()} == dict(widths)
            assert {int(k):n for k,n in expected['raw_START_minus_previous_END_histogram'].items()} == dict(deltas)
            assert expected['extent_checks'] == dict(zero_based_half_open_valid=valid0,one_based_closed_valid=valid1)
            assert expected['contiguous_source_blocks'] == 1
            closed_adjacent = rows[0][1] == 1 and (not deltas or set(deltas) == {1})
            assert closed_adjacent
            covered_end = rows[-1][2]
            if covered_end < lengths[name]:
                omitted.append(dict(chromosome=name,start0=covered_end,end0=lengths[name],
                    status='Outside published table; not certified neutral'))
            summary[chrom] = dict(intervals=len(rows),first_start=rows[0][1],last_end=covered_end,
                FAI_length=lengths[name],one_based_closed_adjacent=True,full_contig_extent=covered_end==lengths[name])
        for name in [str(i) for i in range(1,23)]+['X','Y','MT']:
            if 'chr'+('M' if name=='MT' else name) not in per_chrom and name not in per_chrom:
                omitted.append(dict(chromosome=name,start0=0,end0=lengths[name],
                    status='Chromosome absent from published table; not certified neutral'))
        with (args.out/'intervals_zero_based.tsv').open('x') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['source_row','chromosome','start0','end0','feature'])
            writer.writerows(zero_based)
        (args.out/'report.json').write_text(json.dumps(dict(all_coordinate_rows_exact=len(observed),
            primary_interval_source_files_verified=len(sources),chromosomes=summary,
            source_coordinate_interpretation='One-based closed; source-backed interpretation consistent with exact adjacency',
            converted_interval_unit='Zero-based half-open, start=START-1 and end=END',
            physical_width_histogram=dict(sorted(width_counts.items())),outside_published_table=omitted,
            fixed_500kb_intervals=width_counts[500000],
            original_reference_sequence_identity='Not established',
            current_Ginkgo_code_proves_executed_published_version=False,
            unique_mappability_or_callability_mask='Not established by interval extent',
            assay_values_parsed=False, DNA_CNA_truth_constructed=False,
            biological_null_specimens_verified=0), indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__': main()
