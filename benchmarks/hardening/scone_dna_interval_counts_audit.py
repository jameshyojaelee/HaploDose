"""Independent closed-coordinate reconstruction of published-interval DNA support."""
import argparse
from collections import Counter, defaultdict
import csv
import json
import os
from pathlib import Path
import subprocess

import numpy as np

from haplodose.provenance import start_run, finish_run, file_record


def assign(chromosome, position1, regions):
    if chromosome not in regions:
        return '@outside:chromosome_absent'
    rows = regions[chromosome]
    index = int(np.searchsorted(rows['ends'], position1, side='left'))
    if index == len(rows['ends']):
        return '@outside:after_table'
    if position1 < rows['starts'][index]:
        return '@outside:internal_gap'
    return rows['features'][index]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    chain = args.task.parent/'scone_count_chain_20261109'
    inputs = dict(script=Path(__file__), design=args.task/'INTERVAL_COUNT_AUDIT_DESIGN.md',
        counts=args.task/'interval_counts_v1/run_manifest.json',
        coordinates=args.task/'coordinates_v1/run_manifest.json',
        tools=chain/'tool_environment_v1/run_manifest.json',
        mapping=chain/'dna_map_v1/run_manifest.json')
    run = start_run(args.out, inputs, {}, dict(scope='Independent integer and identity audit; no CNA truth'),
                    mode='scone_dna_interval_counts_audit')
    try:
        for name, path in inputs.items():
            if name in ('script', 'design'): continue
            prior = json.loads(path.read_text()); assert prior['status'] == 'complete'
            for record in prior['output']['files'].values():
                assert file_record(record['path'])['sha256'] == record['sha256']
        fixture = {'1': dict(starts=np.array([1, 11]), ends=np.array([10, 20]), features=['a', 'b'])}
        assert [assign('1', p, fixture) for p in (1, 10, 11, 20, 21)] == ['a', 'a', 'b', 'b', '@outside:after_table']
        assert assign('Y', 1, fixture) == '@outside:chromosome_absent'
        assert assign('1', 5, {'1': dict(starts=np.array([11]), ends=np.array([20]), features=['b'])}) == '@outside:internal_gap'
        original = list(csv.DictReader((args.task/'coordinates_v1/coordinates.tsv').open(), delimiter='\t'))
        grouped = defaultdict(list)
        for row in original:
            chrom = row['CHR_literal'].removeprefix('chr')
            grouped[chrom].append(row)
        regions = {chrom: dict(starts=np.array([int(r['START_literal']) for r in rows]),
            ends=np.array([int(r['END_literal']) for r in rows]), features=[r['feature_literal'] for r in rows])
            for chrom, rows in grouped.items()}
        final_x = regions['X']['features'][-1]
        output = args.task/'interval_counts_v1'
        tools = json.loads((chain/'tool_environment_v1/report.json').read_text()); prefix = Path(tools['environment'])
        samtools = prefix/'bin/samtools'
        assert file_record(samtools)['sha256'] == tools['versions']['samtools']['executable']['sha256']
        environment = {**os.environ, 'LD_LIBRARY_PATH': str(prefix/'lib')}
        mates = defaultdict(dict); starts = Counter(); strata = Counter(); outside = Counter(); states = Counter()
        panels = defaultdict(list); x_quality = Counter()
        with (output/'read_targets.tsv').open() as table, (args.out/'samtools.stderr').open('w') as error:
            exported = csv.DictReader(table, delimiter='\t')
            process = subprocess.Popen([str(samtools), 'view', '-F', '2304', str(chain/'dna_map_v1/SRR16194948/DNA.bam')],
                env=environment, stdout=subprocess.PIPE, stderr=error, text=True)
            try:
                for line in process.stdout:
                    fields = line.rstrip('\n').split('\t'); flag = int(fields[1]); pos = int(fields[3]); quality = int(fields[4])
                    assert bool(flag & 64) != bool(flag & 128)
                    mate = 1 if flag & 64 else 2; name = fields[0]; chrom = fields[2]
                    assert mate not in mates[name]
                    target = '@unmapped' if flag & 4 else assign(chrom, pos, regions)
                    mates[name][mate] = (flag, chrom, pos, quality, fields[5], target)
                    row = next(exported)
                    assert row == dict(QNAME=name, mate=str(mate), chromosome=chrom, POS1=str(pos),
                        MAPQ=str(quality), flag=str(flag), CIGAR=fields[5], target=target)
                    if flag & 4:
                        states['unmapped'] += 1
                        continue
                    states['mapped'] += 1; starts[target] += 1; strata['mapped_mate', target, quality] += 1
                    if target.startswith('@outside:'): outside[chrom, target] += 1
                    if target == '@outside:after_table': assert chrom == 'X'
                    label = None
                    if chrom == 'X' and target == '@outside:after_table':
                        label = 'X_after_table'; x_quality[quality] += 1
                    elif target == final_x: label = 'final_included_X_bin'
                    if label:
                        panels[label].append(dict(QNAME=name, mate=mate, flag=flag, chromosome=chrom, POS1=pos,
                            MAPQ=quality, CIGAR=fields[5], sequence=fields[9], quality=fields[10], target=target))
                        panels[label].sort(key=lambda r: (r['QNAME'], r['mate'])); del panels[label][10:]
                assert next(exported, None) is None
                assert process.wait() == 0
            finally:
                process.stdout.close()
                if process.poll() is None: process.terminate(); process.wait()
        proper = Counter(); candidates = defaultdict(set); crossing = 0
        with (output/'template_targets.tsv').open() as table:
            exported = csv.DictReader(table, delimiter='\t')
            for name in sorted(mates):
                assert set(mates[name]) == {1, 2}
                pair = [mates[name][1], mates[name][2]]
                eligible = not any(r[0] & 4 for r in pair) and all(r[0] & 2 for r in pair) and pair[0][1] == pair[1][1]
                left = min(r[2] for r in pair); quality = min(r[3] for r in pair)
                target = assign(pair[0][1], left, regions) if eligible else '@excluded:nonproper_or_unmapped'
                fingerprint = (tuple((r[1], r[2], bool(r[0] & 16), r[4]) for r in pair), name.rsplit(':umi_', 1)[1])
                assert next(exported) == dict(QNAME=name, proper_same_chr=str(int(eligible)), left_POS1=str(left),
                    min_MAPQ=str(quality), target=target, coordinate_UMI_fingerprint=json.dumps(fingerprint, separators=(',', ':')))
                if eligible:
                    proper[target] += 1; candidates[target].add(fingerprint)
                    strata['proper_template', target, quality] += 1
                    crossing += pair[0][5] != pair[1][5]
            assert next(exported, None) is None
        with (output/'interval_counts.tsv').open() as table:
            exported = csv.DictReader(table, delimiter='\t')
            for source in original:
                feature = source['feature_literal']; row = next(exported)
                assert row == dict(source_row=source['source_row'], feature=feature,
                    chromosome=source['CHR_literal'].removeprefix('chr'), start0=str(int(source['START_literal'])-1),
                    end0=source['END_literal'], mapped_mate_starts=str(starts[feature]),
                    proper_template_left_starts=str(proper[feature]), coordinate_UMI_fingerprints=str(len(candidates[feature])))
            assert next(exported, None) is None
        exported_strata = {}
        with (output/'MAPQ_strata.tsv').open() as table:
            for row in csv.DictReader(table, delimiter='\t'):
                key = (row['unit'], row['target'], int(row['MAPQ']))
                assert key not in exported_strata; exported_strata[key] = int(row['count'])
        assert exported_strata == dict(strata)
        assert json.loads((output/'read_panels.json').read_text()) == dict(panels)
        report = json.loads((output/'report.json').read_text())
        assert report['primary_templates'] == len(mates)
        assert report['primary_mate_states'] == dict(states)
        assert report['mapped_mate_starts'] == sum(starts.values())
        assert report['mapped_starts_inside_published'] == sum(n for k, n in starts.items() if not k.startswith('@'))
        assert report['proper_template_left_starts'] == sum(proper.values())
        assert report['proper_starts_inside_published'] == sum(n for k, n in proper.items() if not k.startswith('@'))
        assert report['coordinate_UMI_fingerprints'] == len(set().union(*candidates.values()))
        assert report['proper_pairs_with_different_start_targets'] == crossing
        assert report['outside_mapped_mates'] == [dict(chromosome=c, reason=t, count=n) for (c, t), n in sorted(outside.items())]
        assert report['outside_proper_template_starts'] == {k:n for k,n in proper.items() if k.startswith('@')}
        assert report['X_after_table_MAPQ_histogram'] == {str(k):v for k,v in x_quality.items()}
        result = dict(fixtures_passed=7, all_primary_mates_exact=sum(len(v) for v in mates.values()),
            all_templates_exact=len(mates), all_interval_rows_exact=len(original), all_MAPQ_strata_exact=len(strata),
            all_panel_records_exact=sum(map(len, panels.values())), all_report_counts_exact=True,
            X_after_table_mapped_mates=sum(x_quality.values()), X_after_table_MAPQ60=x_quality[60],
            original_deduplication_reproduced=False, DNA_CNA_truth_constructed=False, biological_null_specimens_verified=0)
        (args.out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__': main()
