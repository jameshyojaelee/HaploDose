"""Independent DNA primary-pair, original sequence/quality and assembly audit."""
import argparse
from collections import Counter
import gzip
import json
import os
from pathlib import Path
import subprocess

from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    inputs=dict(script=Path(__file__),design=args.task/'DNA_MAP_DESIGN.md')
    for name in ('dna_map_v1','prepare_v1','dna_index_v2','tool_environment_v1'):
        inputs[name]=args.task/name/'run_manifest.json'
    run=start_run(args.out,inputs,{},dict(scope='Exact molecular record/assembly correspondence, not DNA biological truth'),mode='scone_dna_map_audit')
    try:
        for name,path in inputs.items():
            if name in ('script','design'):continue
            manifest=json.loads(path.read_text());assert manifest['status']=='complete'
            for item in manifest['output']['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
        mapped=json.loads((args.task/'dna_map_v1/report.json').read_text());prepared=json.loads((args.task/'prepare_v1/report.json').read_text())
        index=json.loads((args.task/'dna_index_v2/report.json').read_text());fai=Path(index['private_reference']['path']+'.fai')
        expected_sq={line.split('\t')[0]:int(line.split('\t')[1]) for line in fai.read_text().splitlines()}
        tool=json.loads((args.task/'tool_environment_v1/report.json').read_text());prefix=Path(tool['environment'])
        assert file_record(prefix/'bin/samtools')['sha256']==tool['versions']['samtools']['executable']['sha256']
        env={**os.environ,'LD_LIBRARY_PATH':str(prefix/'lib')};results={}
        complement=str.maketrans('ACGTNacgtn','TGCANtgcan')
        for accession in ('SRR16194948','SRR22047451'):
            source=prepared['libraries'][accession];expected={}
            for item in source['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
            for mate in (1,2):
                with gzip.open(source['files'][f'DNA_R{mate}.fastq.gz']['path'],'rt',encoding='ascii') as handle:
                    while True:
                        header=handle.readline()
                        if not header:break
                        sequence=handle.readline().rstrip('\r\n');plus=handle.readline();quality=handle.readline().rstrip('\r\n')
                        assert header.startswith('@') and plus.startswith('+') and len(sequence)==len(quality)
                        key=(header[1:].split()[0],mate);assert key not in expected;expected[key]=(sequence,quality)
            assert len(expected)==2*source['counts']['DNA']
            for item in mapped['libraries'][accession]['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
            bam=args.task/'dna_map_v1'/accession/'DNA.bam'
            header=subprocess.run([str(prefix/'bin/samtools'),'view','-H',str(bam)],env=env,capture_output=True,text=True,timeout=30,check=True)
            observed_sq={}
            for line in header.stdout.splitlines():
                if line.startswith('@SQ\t'):
                    fields=dict(field.split(':',1) for field in line.split('\t')[1:]);assert fields['SN'] not in observed_sq
                    observed_sq[fields['SN']]=int(fields['LN'])
            assert observed_sq==expected_sq
            seen=set();states=Counter();mapq=Counter();chromosomes=Counter();templates={}
            with (args.out/(accession+'.samtools.stderr')).open('w') as error:
                process=subprocess.Popen([str(prefix/'bin/samtools'),'view','-F','2304',str(bam)],env=env,
                    stdout=subprocess.PIPE,stderr=error,text=True)
                try:
                    for line in process.stdout:
                        row=line.rstrip('\n').split('\t');flag=int(row[1]);assert bool(flag&64)!=bool(flag&128)
                        key=(row[0],1 if flag&64 else 2);assert key in expected and key not in seen;seen.add(key)
                        sequence,quality=row[9:11]
                        if flag&16:sequence=sequence.translate(complement)[::-1];quality=quality[::-1]
                        assert (sequence,quality)==expected[key],f'Primary alignment lost original sequence/quality: {key}'
                        state='unmapped' if flag&4 else 'mapped';states[state]+=1
                        if state=='mapped':mapq[int(row[4])]+=1;chromosomes[row[2]]+=1
                        templates.setdefault(row[0],[]).append((flag,row[2],int(row[4])))
                    assert process.wait()==0
                finally:
                    process.stdout.close()
                    if process.poll() is None:process.terminate();process.wait()
            assert seen==set(expected) and all(len(value)==2 for value in templates.values())
            paired=Counter();pair_q=Counter()
            for rows in templates.values():
                if any(flag&4 for flag,chrom,q in rows):paired['one_or_both_unmapped']+=1
                elif rows[0][1]!=rows[1][1]:paired['different_chromosomes']+=1
                elif all(flag&2 for flag,chrom,q in rows):
                    paired['both_mapped_same_chr_proper']+=1;pair_q[min(row[2] for row in rows)]+=1
                else:paired['both_mapped_same_chr_nonproper']+=1
            assert sum(paired.values())==source['counts']['DNA'] and sum(mapq.values())==states['mapped']
            results[accession]=dict(primary_templates=len(templates),primary_records=len(seen),states=dict(states),
                original_mate_sequences_qualities_exact=True,reference_contigs_exact=len(observed_sq),
                mapped_read_MAPQ_histogram=dict(sorted(mapq.items())),mapped_reads_by_chromosome=dict(chromosomes),
                primary_pair_states=dict(paired),proper_pair_min_MAPQ_histogram=dict(sorted(pair_q.items())))
        (args.out/'report.json').write_text(json.dumps(dict(libraries=results,MAPQ_filter=None,PCR_UMI_deduplicated=False,
            DNA_copy_number_truth=False,biological_null_specimens_verified=0),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
