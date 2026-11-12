"""Independent SAM-to-table audit of the entire locked DNA support census."""
import argparse
from collections import Counter, defaultdict
import csv
import json
import os
from pathlib import Path
import subprocess

from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'DNA_SUPPORT_DESIGN.md')
    for name in ('dna_support_v1', 'dna_map_v1', 'dna_map_audit_v1', 'dna_index_v2', 'tool_environment_v1'):
        inputs[name] = args.task/name/'run_manifest.json'
    run = start_run(args.out, inputs, {}, dict(scope='Independent exact integer/identity reference; no DNA truth'),
                    mode='scone_dna_support_audit')
    try:
        for name, path in inputs.items():
            if name in ('script', 'design'): continue
            manifest = json.loads(path.read_text()); assert manifest['status'] == 'complete'
            for record in manifest['output']['files'].values():
                assert file_record(record['path'])['sha256'] == record['sha256']
        report = json.loads((args.task/'dna_support_v1/report.json').read_text())
        tools = json.loads((args.task/'tool_environment_v1/report.json').read_text())
        prefix = Path(tools['environment']); env = {**os.environ, 'LD_LIBRARY_PATH':str(prefix/'lib')}
        assert file_record(prefix/'bin/samtools')['sha256'] == tools['versions']['samtools']['executable']['sha256']
        index = json.loads((args.task/'dna_index_v2/report.json').read_text())
        fai = Path(index['private_reference']['path']+'.fai')
        lengths = {line.split('\t')[0]:int(line.split('\t')[1]) for line in fai.read_text().splitlines()}
        exported = defaultdict(dict)
        with (args.task/'dna_support_v1/template_fingerprints.tsv').open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                assert row['QNAME'] not in exported[row['library']]
                exported[row['library']][row['QNAME']] = row
        assert set(exported) == {'SRR16194948', 'SRR22047451'}
        results = {}
        for accession, expected in exported.items():
            primary = defaultdict(dict)
            bam = args.task/'dna_map_v1'/accession/'DNA.bam'
            with (args.out/(accession+'.samtools.stderr')).open('w') as error:
                process = subprocess.Popen([str(prefix/'bin/samtools'), 'view', '-F', '2304', str(bam)],
                    env=env, stdout=subprocess.PIPE, stderr=error, text=True)
                try:
                    for line in process.stdout:
                        row = line.rstrip('\n').split('\t'); flag = int(row[1])
                        mate = 1 if flag & 0x40 else 2
                        assert bool(flag&0x40) != bool(flag&0x80)
                        assert row[0] in expected and mate not in primary[row[0]]
                        primary[row[0]][mate] = (flag, row[2], int(row[3]), int(row[4]), row[5])
                    assert process.wait() == 0
                finally:
                    process.stdout.close()
                    if process.poll() is None: process.terminate(); process.wait()
            assert set(primary) == set(expected)
            coordinates = Counter(); molecular_candidates = Counter(); strata = defaultdict(Counter)
            bin_counts = defaultdict(Counter); umi_lengths = Counter(); sequence_types = Counter(); total = 0
            for name, mates in primary.items():
                assert set(mates) == {1,2}
                first, second = mates[1], mates[2]
                key = tuple((value[1],value[2],bool(value[0]&0x10),value[4]) for value in (first,second))
                encoded = json.dumps(key, separators=(',', ':'))
                source = expected[name]
                assert source['coordinate_fingerprint'] == encoded
                pieces = name.rsplit(':umi_', 1); assert len(pieces) == 2
                umi = pieces[1]; assert source['candidate_UMI'] == umi
                umi_lengths[len(umi)] += 1
                sequence_types['contains_N' if 'N' in umi else 'ACGT_only'] += 1
                quality = min(first[3],second[3]); assert int(source['min_MAPQ']) == quality
                eligible = not(first[0]&0x4 or second[0]&0x4) and bool(first[0]&0x2 and second[0]&0x2) and first[1] == second[1]
                assert bool(int(source['primary_proper_same_chr'])) == eligible
                if not eligible: continue
                total += 1; coordinates[encoded] += 1; molecular_candidates[encoded,umi] += 1
                strata[quality][encoded,umi] += 1
                for width in (500000,1000000,5000000):
                    # Convert each one-based mate position independently before
                    # choosing the smallest observed coordinate.
                    start = min(first[2]-1, second[2]-1)
                    bin_counts[width,first[1],start//width][encoded,umi] += 1
            observed = report['libraries'][accession]
            assert observed['primary_templates'] == len(primary)
            assert observed['proper_same_chr_templates'] == total
            assert observed['excluded_from_fingerprint_census'] == len(primary)-total
            assert {int(k):v for k,v in observed['candidate_UMI_lengths'].items()} == dict(umi_lengths)
            assert observed['candidate_UMI_sequence_types'] == dict(sequence_types)
            assert observed['coordinate_fingerprints'] == len(coordinates)
            assert observed['coordinate_UMI_fingerprints'] == len(molecular_candidates)
            assert {int(k):v for k,v in observed['coordinate_multiplicity_histogram'].items()} == dict(Counter(coordinates.values()))
            assert {int(k):v for k,v in observed['coordinate_UMI_multiplicity_histogram'].items()} == dict(Counter(molecular_candidates.values()))
            expected_strata = {quality:dict(templates=sum(data.values()),coordinate_fingerprints=len({key for key,umi in data}),
                coordinate_UMI_fingerprints=len(data)) for quality,data in strata.items()}
            assert {int(k):v for k,v in observed['MAPQ_strata'].items()} == expected_strata
            universe = {(width, chrom, number) for width in (500000,1000000,5000000)
                for chrom,length in lengths.items() for number in range((length+width-1)//width)}
            seen = set(); totals = Counter()
            with (args.task/'dna_support_v1'/(accession+'.bins.tsv')).open() as handle:
                for row in csv.DictReader(handle, delimiter='\t'):
                    width = int(row['width']); start = int(row['start0']); chrom = row['chromosome']
                    assert start % width == 0
                    key = (width,chrom,start//width)
                    assert key in universe and key not in seen; seen.add(key)
                    assert int(row['end0']) == min(start+width,lengths[chrom])
                    data = bin_counts.get(key, {})
                    assert int(row['templates']) == sum(data.values())
                    assert int(row['coordinate_fingerprints']) == len({coord for coord,umi in data})
                    assert int(row['coordinate_UMI_fingerprints']) == len(data)
                    totals[width] += int(row['templates'])
            assert seen == universe and set(totals.values()) == {total}
            results[accession] = dict(all_templates_exact=len(primary), all_bin_rows_exact=len(seen),
                all_MAPQ_strata_exact=len(strata), all_multiplicity_histograms_exact=True,
                all_resolution_totals_exact=True)
        (args.out/'report.json').write_text(json.dumps(dict(libraries=results,
            independent_parser=True, CNA_truth_constructed=False, biological_null_specimens_verified=0), indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__': main()
