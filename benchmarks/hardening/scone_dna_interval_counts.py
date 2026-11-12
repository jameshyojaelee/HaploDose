"""Observed DNA start support on exact published intervals, retaining all outside support."""
import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
import csv
import json
import os
from pathlib import Path
import subprocess

from haplodose.provenance import start_run, finish_run, file_record


def locate(chrom, position0, intervals):
    if chrom not in intervals: return '@outside:chromosome_absent'
    rows, ends = intervals[chrom]
    number = bisect_right(ends, position0)
    if number == len(rows): return '@outside:after_table'
    first,last,feature = rows[number]
    if position0 < first: return '@outside:internal_gap'
    assert first <= position0 < last
    return feature


def fixtures():
    intervals = {'1': ([(0,10,'a'),(10,20,'b')],[10,20])}
    expected = ['a','a','b','b','@outside:after_table']
    assert [locate('1',p,intervals) for p in (0,9,10,19,20)] == expected
    assert locate('Y',0,intervals) == '@outside:chromosome_absent'
    assert locate('1',5,{'1': ([(10,20,'b')],[20])}) == '@outside:internal_gap'
    return 7


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True); parser.add_argument('--task',type=Path,required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    counts = args.task.parent/'scone_count_chain_20261109'
    inputs = dict(script=Path(__file__),design=args.task/'INTERVAL_COUNT_DESIGN.md',
        coordinates=args.task/'coordinate_audit_v2/run_manifest.json')
    for name in ('dna_map_v1','dna_map_audit_v1','dna_index_v2','dna_support_v1','tool_environment_v1'):
        inputs[name] = counts/name/'run_manifest.json'
    run = start_run(args.out,inputs,dict(locked_run='SRR16194948',MAPQ_cutoff=None),
        dict(scope='Exact physical support of observed starts, not DNA dosage or original normalized-count reproduction'),
        mode='scone_dna_interval_counts')
    try:
        checked = fixtures()
        for name,path in inputs.items():
            if name in ('script','design'): continue
            prior = json.loads(path.read_text()); assert prior['status'] == 'complete'
            for record in prior['output']['files'].values():
                assert file_record(record['path'])['sha256'] == record['sha256']
        intervals = defaultdict(list); annotation = {}
        with (args.task/'coordinate_audit_v2/intervals_zero_based.tsv').open() as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                first,last = int(row['start0']),int(row['end0']); feature = row['feature']
                assert feature not in annotation
                intervals[row['chromosome']].append((first,last,feature)); annotation[feature] = row
        for chrom,rows in list(intervals.items()):
            assert rows == sorted(rows) and all(a[1] == b[0] for a,b in zip(rows,rows[1:]))
            intervals[chrom] = (rows,[row[1] for row in rows])
        tools = json.loads((counts/'tool_environment_v1/report.json').read_text()); prefix = Path(tools['environment'])
        assert file_record(prefix/'bin/samtools')['sha256'] == tools['versions']['samtools']['executable']['sha256']
        env = {**os.environ,'LD_LIBRARY_PATH':str(prefix/'lib')}
        earlier = json.loads((counts/'dna_map_audit_v1/report.json').read_text())['libraries']['SRR16194948']
        support = json.loads((counts/'dna_support_v1/report.json').read_text())['libraries']['SRR16194948']
        templates = defaultdict(list); starts = Counter(); strata = Counter(); by_chrom = Counter(); outside = Counter()
        panels = defaultdict(list); mate_states = Counter()
        final_X_feature = intervals['X'][0][-1][2]
        bam = counts/'dna_map_v1/SRR16194948/DNA.bam'
        with (args.out/'read_targets.tsv').open('x') as table,(args.out/'samtools.stderr').open('w') as error:
            writer = csv.writer(table,delimiter='\t'); writer.writerow(['QNAME','mate','chromosome','POS1','MAPQ','flag','CIGAR','target'])
            process = subprocess.Popen([str(prefix/'bin/samtools'),'view','-F','2304',str(bam)],
                env=env,stdout=subprocess.PIPE,stderr=error,text=True)
            try:
                for line in process.stdout:
                    row = line.rstrip('\n').split('\t'); flag = int(row[1]); chrom = row[2]; pos = int(row[3]); q = int(row[4])
                    assert bool(flag&64) != bool(flag&128)
                    mate = 1 if flag&64 else 2
                    target = '@unmapped' if flag&4 else locate(chrom,pos-1,intervals)
                    item = dict(name=row[0],mate=mate,chrom=chrom,pos=pos,q=q,flag=flag,cigar=row[5],target=target)
                    templates[row[0]].append(item)
                    writer.writerow([row[0],mate,chrom,pos,q,flag,row[5],target])
                    if flag&4: mate_states['unmapped'] += 1; continue
                    assert pos > 0
                    mate_states['mapped'] += 1; starts[target] += 1; strata['mapped_mate',target,q] += 1
                    by_chrom[chrom] += 1
                    if target.startswith('@outside:'): outside[chrom,target] += 1
                    label = 'X_after_table' if chrom=='X' and target=='@outside:after_table' else 'final_included_X_bin' if target==final_X_feature else None
                    if label:
                        panel = panels[label]
                        panel.append(dict(QNAME=row[0],mate=mate,flag=flag,chromosome=chrom,POS1=pos,MAPQ=q,
                            CIGAR=row[5],sequence=row[9],quality=row[10],target=target))
                        panel.sort(key=lambda value:(value['QNAME'],value['mate']))
                        del panel[10:]
                assert process.wait() == 0
            finally:
                process.stdout.close()
                if process.poll() is None: process.terminate(); process.wait()
        assert len(templates) == earlier['primary_templates']
        assert dict(mate_states) == earlier['states'] and dict(by_chrom) == earlier['mapped_reads_by_chromosome']
        proper = Counter(); candidates = defaultdict(set); global_candidates = {}; crossing = 0; total_proper = 0
        with (args.out/'template_targets.tsv').open('x') as handle:
            writer = csv.writer(handle,delimiter='\t')
            writer.writerow(['QNAME','proper_same_chr','left_POS1','min_MAPQ','target','coordinate_UMI_fingerprint'])
            for name,rows in sorted(templates.items()):
                assert len(rows)==2 and {row['mate'] for row in rows}=={1,2}
                rows.sort(key=lambda row:row['mate'])
                eligible = all(not row['flag']&4 and row['flag']&2 for row in rows) and rows[0]['chrom']==rows[1]['chrom']
                pos = min(row['pos'] for row in rows); q = min(row['q'] for row in rows)
                target = locate(rows[0]['chrom'],pos-1,intervals) if eligible else '@excluded:nonproper_or_unmapped'
                key = tuple((row['chrom'],row['pos'],bool(row['flag']&16),row['cigar']) for row in rows)
                fingerprint = (key,name.rsplit(':umi_',1)[1])
                writer.writerow([name,int(eligible),pos,q,target,json.dumps(fingerprint,separators=(',',':'))])
                if not eligible: continue
                total_proper += 1; proper[target] += 1; candidates[target].add(fingerprint)
                assert fingerprint not in global_candidates or global_candidates[fingerprint]==target
                global_candidates[fingerprint] = target; strata['proper_template',target,q] += 1
                crossing += int(rows[0]['target'] != rows[1]['target'])
        assert total_proper == earlier['primary_pair_states']['both_mapped_same_chr_proper']
        assert len(global_candidates) == support['coordinate_UMI_fingerprints']
        assert sum(starts.values()) == earlier['states']['mapped']
        assert sum(proper.values()) == total_proper
        with (args.out/'interval_counts.tsv').open('x') as handle:
            writer = csv.writer(handle,delimiter='\t')
            writer.writerow(['source_row','feature','chromosome','start0','end0','mapped_mate_starts','proper_template_left_starts','coordinate_UMI_fingerprints'])
            for feature,row in annotation.items():
                writer.writerow([row['source_row'],feature,row['chromosome'],row['start0'],row['end0'],starts[feature],proper[feature],len(candidates[feature])])
        with (args.out/'MAPQ_strata.tsv').open('x') as handle:
            writer = csv.writer(handle,delimiter='\t'); writer.writerow(['unit','target','MAPQ','count'])
            writer.writerows((*key,n) for key,n in sorted(strata.items()))
        summary = dict(fixtures_passed=checked,primary_templates=len(templates),primary_mate_states=dict(mate_states),
            mapped_mate_starts=sum(starts.values()),mapped_starts_inside_published=sum(n for key,n in starts.items() if not key.startswith('@')),
            proper_template_left_starts=total_proper,proper_starts_inside_published=sum(n for key,n in proper.items() if not key.startswith('@')),
            coordinate_UMI_fingerprints=len(global_candidates),proper_pairs_with_different_start_targets=crossing,
            outside_mapped_mates=[dict(chromosome=key[0],reason=key[1],count=n) for key,n in sorted(outside.items())],
            outside_proper_template_starts={key:n for key,n in proper.items() if key.startswith('@')},
            X_after_table_MAPQ_histogram={str(q):strata['mapped_mate','@outside:after_table',q] for q in range(61) if strata['mapped_mate','@outside:after_table',q]},
            observation_unit='Aligned primary mate start or leftmost proper primary pair start, explicitly distinct',
            MAPQ_cutoff=None,original_deduplication_reproduced=False,independent_molecules_certified=False,
            published_assay_values_parsed=False,DNA_CNA_truth_constructed=False,biological_null_specimens_verified=0,
            independent_audit='Not run')
        (args.out/'read_panels.json').write_text(json.dumps(dict(panels),indent=2)+'\n')
        (args.out/'report.json').write_text(json.dumps(summary,indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error)
        raise


if __name__=='__main__': main()
