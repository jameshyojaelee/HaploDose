"""Executed DNA-bin counterexamples and bounded primary SAM evidence."""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import subprocess

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def awk_program(path):
    match=re.search(r" '(?P<body>(?:\{|/\^@RG).*?)' \\\n",path.read_text(),re.S)
    if match is None:
        raise ValueError('Cannot identify executed awk program: '+str(path))
    return match.group('body')


def execute(program,sam,header):
    result=subprocess.run(['awk','-F','\t','-v','OFS=\t','-v','H='+str(header),program],
                          input=sam,text=True,capture_output=True,check=True)
    values=Counter()
    for row in result.stdout.splitlines():
        rg,chrom,bin_id,count=row.split('\t');values[(rg,chrom,int(bin_id))]+=int(count)
    return values,result


def expected_counts(sam):
    counts=Counter();moved=[]
    for row in sam.splitlines():
        if row.startswith('@'):
            continue
        fields=row.split('\t');rg=next((v[5:] for v in fields[11:] if v.startswith('RG:Z:')),None)
        if rg is None:
            continue
        pos=int(fields[3]);zero=pos-1
        # Independent interval-membership reference instead of division.
        candidates=[i for i in range(500) if i*500000<=zero<(i+1)*500000]
        assert len(candidates)==1
        counts[(rg,fields[2],candidates[0])]+=1
        if pos in range(500000,250000001,500000):
            moved.append(dict(read=fields[0],rg=rg,chromosome=fields[2],SAM_POS=pos,
                              zero_based_position=zero,correct_bin=candidates[0]))
    return counts,moved


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--fixtures-only',action='store_true')
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    source=Path(__file__).resolve().parents[2]
    scripts=['fetch_bins.sbatch','fetch_original_bam.sbatch']
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md')
    for name in scripts:
        inputs['old:'+name]=args.task/'baseline/benchmarks/andor'/name
        inputs['new:'+name]=source/'benchmarks/andor'/name
    run=start_run(args.out,inputs,dict(boundaries=0 if args.fixtures_only else 20,
                  fixtures_only=args.fixtures_only,query_timeout_seconds=60),
                  dict(scope='Boundary computation plus bounded reads; historical DNA truth unchanged'),mode='andor_dna_coordinates')
    try:
        versions={}
        for command in (['sam-dump','--version'],['samtools','--version']):
            result=subprocess.run(command,capture_output=True,text=True,check=True)
            versions[command[0]]=result.stdout+result.stderr
        write_json(versions,args.out/'versions.json')
        positions=[1,499999,500000,500001,999999,1000000,1000001]
        sam='@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n@SQ\tSN:chr2\tLN:242193529\n@RG\tID:A\tSM:cellA\n@RG\tID:B\tSM:cellB\n'
        for i,pos in enumerate(positions):
            sam+=f'r{i}\t99\tchr1\t{pos}\t60\t1M\t=\t{pos+2}\t3\tA\tI\tRG:Z:A\n'
        sam+='other\t99\tchr2\t500000\t60\t1M\t=\t500002\t3\tA\tI\tRG:Z:B\n'
        sam+='missing\t99\tchr1\t500000\t60\t1M\t=\t500002\t3\tA\tI\n'
        (args.out/'fixture.sam').write_text(sam)
        expected,moved=expected_counts(sam);checks=[]
        assert len(moved)==3 and sum(expected.values())==8
        programs={}
        for name in scripts:
            outputs={}
            for version in ('old','new'):
                program=awk_program(inputs[version+':'+name]);programs[(version,name)]=program
                result,process=execute(program,sam,args.out/(version+'_'+name+'.headers'))
                outputs[version]=result
                (args.out/(version+'_'+name+'.stdout')).write_text(process.stdout)
                (args.out/(version+'_'+name+'.stderr')).write_text(process.stderr)
            assert outputs['new']==expected and outputs['old']!=expected
            assert sum(outputs['new'].values())==sum(outputs['old'].values())==8
            changed=0
            headers='\n'.join(row for row in sam.splitlines() if row.startswith('@'))+'\n'
            for row in sam.splitlines():
                if row.startswith('@'):
                    continue
                single=headers+row+'\n';answer,_=expected_counts(single)
                before,_=execute(programs[('old',name)],single,args.out/'single_old.headers')
                after,_=execute(programs[('new',name)],single,args.out/'single_new.headers')
                assert after==answer and sum(before.values())==sum(after.values())
                changed+=before!=after
            assert changed==3
            checks.append(dict(script=name,filtered_fragments=8,corrected_boundary_assignments=changed,
                               net_positive_bin_count_change=sum((outputs['new']-outputs['old']).values()),
                               exact_reference=True))
        write_json(checks,args.out/'fixture_checks.json')
        if args.fixtures_only:
            write_json(dict(fixtures=checks,archive_queries=0,historical_truth_changed=False),args.out/'report.json')
            finish_run(args.out,run,backend='cpu')
            return
        data=Path('/gpfs/commons/home/jameslee/Haplodose/data/external/andor2020_scdna/sra')
        queries=[]
        for line,accession in [('SNU-601','SRR8117283'),('MKN-45','SRR8117281'),('NUGC-4','SRR8117282')]:
            archive=data/accession/(accession+'.sra');stat=archive.stat()
            folder=args.out/line;folder.mkdir()
            for i in range(1,21):
                boundary=i*500000;region=f'chr1:{boundary-1}-{boundary+1}'
                command=['sam-dump','--primary','--aligned-region',region,str(archive)]
                row=dict(line=line,accession=accession,region=region,command=command,
                         archive=dict(path=str(archive),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns,
                                      integrity='full archive checksum not reassessed in bounded inspection'))
                try:
                    result=subprocess.run(command,capture_output=True,text=True,timeout=60)
                    raw=folder/f'{i:02d}.sam';raw.write_text(result.stdout)
                    error=folder/f'{i:02d}.stderr';error.write_text(result.stderr)
                    row.update(exit_code=result.returncode,raw=file_record(raw),stderr=file_record(error))
                    if result.returncode==0:
                        filtered=subprocess.run(['samtools','view','-h','-f','0x42','-F','0xF04','-q','30','-'],
                                                input=result.stdout,text=True,capture_output=True,check=True)
                        path=folder/f'{i:02d}.filtered.sam';path.write_text(filtered.stdout)
                        expected,observed_moves=expected_counts(filtered.stdout)
                        actual,_=execute(programs[('new',scripts[0])],filtered.stdout,folder/f'{i:02d}.headers')
                        assert expected==actual
                        row.update(filtered=file_record(path),filtered_fragments=sum(expected.values()),
                                   observed_boundary_starts=observed_moves)
                    else:
                        row['status']='source_query_failed'
                except subprocess.TimeoutExpired:
                    row.update(status='source_query_timeout',exit_code=None)
                queries.append(row)
                write_json(queries,args.out/'queries.json')
            print(json.dumps(dict(line=line,queries=len([r for r in queries if r['line']==line]),status='complete')),flush=True)
        successful=[r for r in queries if r.get('exit_code')==0]
        moves=[dict(line=r['line'],region=r['region'],**move) for r in successful for move in r['observed_boundary_starts']]
        write_json(dict(fixtures=checks,queries=len(queries),successful_queries=len(successful),
                        source_failures=len(queries)-len(successful),observed_boundary_starts=moves,
                        biological_specimens=3,scope='Development read trace, not new independent validation',
                        historical_truth_changed=False,clone_impact='not established by bounded loci'),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
