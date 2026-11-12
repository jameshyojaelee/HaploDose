"""DNA template/fingerprint support census before constructing any CNA reference."""
import argparse
from collections import Counter, defaultdict
import csv
import json
import os
from pathlib import Path
import re
import subprocess

from haplodose.provenance import start_run, finish_run, file_record


def fingerprint(rows):
    assert len(rows) == 2
    first, second = sorted(rows, key=lambda row: 0 if row['flag'] & 64 else 1)
    assert first['flag'] & 64 and second['flag'] & 128
    # Observed aligned endpoints, orientation and CIGAR only. This is a
    # repeat-record diagnostic, not a molecular identity or duplicate caller.
    return tuple((row['chrom'], row['pos'], bool(row['flag'] & 16), row['cigar'])
                 for row in (first, second))


def fixtures():
    first = dict(flag=99, chrom='1', pos=500000, cigar='10M')
    second = dict(flag=147, chrom='1', pos=500020, cigar='10M')
    assert fingerprint([first, second]) == fingerprint([second, first])
    base = fingerprint([first, second])
    different = fingerprint([first, dict(second, pos=500021)])
    coordinates = Counter([base, base, base, different])
    candidates = Counter([(base, 'AAAAA'), (base, 'AAAAA'), (base, 'AAAAC'), (different, 'AAAAA')])
    assert sum(coordinates.values()) == sum(candidates.values()) == 4
    assert len(coordinates) == 2 and len(candidates) == 3
    assert (500000-1)//500000 == 0 and (500001-1)//500000 == 1
    return 4


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'DNA_SUPPORT_DESIGN.md')
    for name in ('dna_map_v1', 'dna_map_audit_v1', 'dna_index_v2', 'prepare_v1', 'tool_environment_v1'):
        inputs[name] = args.task/name/'run_manifest.json'
    run = start_run(args.out, inputs, dict(bin_widths=[500000,1000000,5000000], MAPQ_cutoff=None),
        dict(scope='Conditional depth/record census; no CNA reference or calibrated molecule count'),
        mode='scone_dna_support')
    try:
        checked = fixtures()
        for name, path in inputs.items():
            if name in ('script', 'design'): continue
            manifest = json.loads(path.read_text()); assert manifest['status'] == 'complete'
            for record in manifest['output']['files'].values():
                assert file_record(record['path'])['sha256'] == record['sha256']
        tools = json.loads((args.task/'tool_environment_v1/report.json').read_text())
        prefix = Path(tools['environment']); env = {**os.environ, 'LD_LIBRARY_PATH': str(prefix/'lib')}
        assert file_record(prefix/'bin/samtools')['sha256'] == tools['versions']['samtools']['executable']['sha256']
        index = json.loads((args.task/'dna_index_v2/report.json').read_text())
        fai = Path(index['private_reference']['path']+'.fai')
        lengths = {row.split('\t')[0]:int(row.split('\t')[1]) for row in fai.read_text().splitlines()}
        earlier = json.loads((args.task/'dna_map_audit_v1/report.json').read_text())
        prepared = json.loads((args.task/'prepare_v1/report.json').read_text())
        reports = {}; all_members = args.out/'template_fingerprints.tsv'
        with all_members.open('x') as output:
            writer = csv.writer(output, delimiter='\t')
            writer.writerow(['library', 'QNAME', 'candidate_UMI', 'min_MAPQ', 'primary_proper_same_chr', 'coordinate_fingerprint'])
            for accession in ('SRR16194948', 'SRR22047451'):
                bam = args.task/'dna_map_v1'/accession/'DNA.bam'
                records = defaultdict(list)
                with (args.out/(accession+'.samtools.stderr')).open('w') as error:
                    process = subprocess.Popen([str(prefix/'bin/samtools'), 'view', '-F', '2304', str(bam)],
                        env=env, stdout=subprocess.PIPE, stderr=error, text=True)
                    try:
                        for line in process.stdout:
                            fields = line.rstrip('\n').split('\t')
                            records[fields[0]].append(dict(flag=int(fields[1]), chrom=fields[2], pos=int(fields[3]),
                                mapq=int(fields[4]), cigar=fields[5]))
                        assert process.wait() == 0
                    finally:
                        process.stdout.close()
                        if process.poll() is None: process.terminate(); process.wait()
                assert len(records) == earlier['libraries'][accession]['primary_templates']
                coordinates = Counter(); candidates = Counter(); rows_by_quality = Counter()
                sequence_types = Counter(); quality_counts = Counter(); binned = defaultdict(Counter)
                proper = 0; excluded = 0; umi_lengths = Counter()
                for name, rows in sorted(records.items()):
                    assert len(rows) == 2 and sum(bool(row['flag']&64) for row in rows) == 1
                    assert sum(bool(row['flag']&128) for row in rows) == 1
                    umi = name.rsplit(':umi_', 1)[1]
                    assert re.fullmatch('[ACGTN]+', umi), 'Unexpected candidate UMI alphabet'
                    umi_lengths[len(umi)] += 1
                    sequence_types['contains_N' if 'N' in umi else 'ACGT_only'] += 1
                    eligible = all(not row['flag']&4 and row['flag']&2 for row in rows) and rows[0]['chrom'] == rows[1]['chrom']
                    q = min(row['mapq'] for row in rows)
                    key = fingerprint(rows)
                    writer.writerow([accession, name, umi, q, int(eligible), json.dumps(key, separators=(',', ':'))])
                    if not eligible: excluded += 1; continue
                    proper += 1; coordinates[key] += 1; candidates[key, umi] += 1
                    rows_by_quality[q, key, umi] += 1; quality_counts[q] += 1
                    chrom = rows[0]['chrom']; left = min(row['pos'] for row in rows)
                    assert 1 <= left <= lengths[chrom]
                    for width in (500000,1000000,5000000):
                        binned[width, chrom, (left-1)//width][q, key, umi] += 1
                assert proper == earlier['libraries'][accession]['primary_pair_states']['both_mapped_same_chr_proper']
                assert proper+excluded == len(records)
                assert sum(coordinates.values()) == sum(candidates.values()) == sum(quality_counts.values()) == proper
                assert len(coordinates) <= len(candidates) <= proper
                # MAPQ is retained as a full categorical diagnostic. A fingerprint
                # seen at several qualities is one fingerprint globally, so these
                # stratum-specific unique counts must not be summed as global n.
                strata = {}
                for q in sorted(quality_counts):
                    data = {(key, umi):n for (quality,key,umi),n in rows_by_quality.items() if quality == q}
                    strata[q] = dict(templates=sum(data.values()), coordinate_fingerprints=len({key for key,umi in data}),
                        coordinate_UMI_fingerprints=len(data))
                bins_path = args.out/(accession+'.bins.tsv')
                with bins_path.open('x') as handle:
                    table = csv.writer(handle, delimiter='\t')
                    table.writerow(['width', 'chromosome', 'start0', 'end0', 'templates', 'coordinate_fingerprints', 'coordinate_UMI_fingerprints'])
                    for width in (500000,1000000,5000000):
                        total = 0
                        for chrom, length in lengths.items():
                            for number in range((length+width-1)//width):
                                data = binned.get((width, chrom, number), {})
                                n = sum(data.values()); total += n
                                table.writerow([width, chrom, number*width, min((number+1)*width,length), n,
                                    len({key for q,key,umi in data}), len({(key,umi) for q,key,umi in data})])
                        assert total == proper
                reports[accession] = dict(primary_templates=len(records), proper_same_chr_templates=proper,
                    excluded_from_fingerprint_census=excluded, candidate_UMI_lengths=dict(umi_lengths),
                    candidate_UMI_sequence_types=dict(sequence_types), coordinate_fingerprints=len(coordinates),
                    coordinate_UMI_fingerprints=len(candidates),
                    coordinate_multiplicity_histogram=dict(Counter(coordinates.values())),
                    coordinate_UMI_multiplicity_histogram=dict(Counter(candidates.values())),
                    MAPQ_strata=strata, all_three_bin_template_totals_exact=True)
        (args.out/'report.json').write_text(json.dumps(dict(fixtures_passed=checked, libraries=reports,
            bin_widths=[500000,1000000,5000000], MAPQ_cutoff=None,
            independent_audit='Not run', UMI_fingerprints_are_certified_molecules=False,
            CNA_truth_constructed=False, biological_null_specimens_verified=0), indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__': main()
