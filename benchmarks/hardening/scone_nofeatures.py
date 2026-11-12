"""Prespecified read-level exon/gene-body mechanism panel, independent of CNA outcomes."""
import argparse
from bisect import bisect_left
from collections import Counter, defaultdict
import csv
import fcntl
import gzip
import json
import os
from pathlib import Path
import re
import subprocess

from haplodose.provenance import start_run, finish_run, file_record


def make_index(intervals):
    ordered=sorted(intervals);ends=[];maximum=-1
    for start,end,gene in ordered:
        assert 0<=start<end
        maximum=max(maximum,end);ends.append(maximum)
    return ordered,[item[0] for item in ordered],ends


def overlap(index,start,end):
    intervals,starts,maximum=index;position=bisect_left(starts,end)-1;genes=set()
    while position>=0 and maximum[position]>start:
        left,right,gene=intervals[position]
        if left<end and right>start:genes.add(gene)
        position-=1
    return genes


def blocks(pos,cigar,query_length):
    pieces=re.findall(r'(\d+)([MIDNSHP=X])',cigar)
    assert ''.join(n+op for n,op in pieces)==cigar
    current=pos-1;query=0;result=[]
    for text,op in pieces:
        length=int(text);assert length>0
        if op in 'M=X':result.append((current,current+length))
        if op in 'MDN=X':current+=length
        if op in 'MIS=X':query+=length
    assert query==query_length
    return result


def count_table(path):
    with path.open() as handle:
        table=list(csv.reader((line for line in handle if not line.startswith('#')),delimiter='\t'))
    assert table[0][:6]==['Geneid','Chr','Start','End','Strand','Length']
    return [row[:6]+[int(row[-1])] for row in table[1:]]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    inputs=dict(script=Path(__file__),design=args.task/'NOFEATURES_DESIGN.md')
    for name in ('rna_map_v1','rna_map_audit_v1','reference_v3','tool_environment_v1'):
        inputs[name]=args.task/name/'run_manifest.json'
    run=start_run(args.out,inputs,dict(panel_sizes={'Unassigned_NoFeatures':50,'Assigned':10,'Unassigned_Ambiguity':10}),
                  dict(scope='Purposive mechanism panel, not population validation'),mode='scone_nofeatures')
    try:
        for name,path in inputs.items():
            if name in ('script','design'):continue
            manifest=json.loads(path.read_text());assert manifest['status']=='complete'
            for item in manifest['output']['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
        assert blocks(101,'2S3M5N2M1I2M1D1M',11)==[(100,103),(108,110),(110,112),(113,114)]
        fixture=[(100,200,'A'),(150,151,'B'),(300,310,'C')];index=make_index(fixture);controls=0
        for start,end in [(99,100),(100,101),(199,200),(200,201),(150,151),(200,300),(0,400)]:
            assert overlap(index,start,end)=={g for a,b,g in fixture if a<end and b>start};controls+=1
        reference=json.loads((args.task/'reference_v3/report.json').read_text())
        tool=json.loads((args.task/'tool_environment_v1/report.json').read_text());prefix=Path(tool['environment'])
        env={**os.environ,'LD_LIBRARY_PATH':str(prefix/'lib')};mapping=json.loads((args.task/'rna_map_v1/report.json').read_text())
        for item in tool['versions'].values():assert file_record(item['executable']['path'])['sha256']==item['executable']['sha256']
        inventory=json.loads((args.task/'tool_environment_v1/tool_files.json').read_text())
        for item in inventory.values():assert file_record(item['path'])['sha256']==item['sha256']
        gtf=Path(reference['reference']['Homo_sapiens.GRCh38.112.gtf.gz']['path']);results={};panels={};commands=[]
        with (gtf.parent.parent/'.setup_reference.lock').open('rb') as lock:
            fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
            assert file_record(gtf)['sha256']==reference['reference'][gtf.name]['sha256']
            raw_annotation=defaultdict(list)
            with gzip.open(gtf,'rt') as handle:
                for line in handle:
                    if line.startswith('#'):continue
                    fields=line.rstrip('\n').split('\t')
                    if fields[2] not in ('exon','gene'):continue
                    match=re.search(r'(?:^|;)\s*gene_id "([^"]+)"',fields[8]);assert match
                    raw_annotation[fields[0],fields[2]].append((int(fields[3])-1,int(fields[4]),match.group(1)))
            annotation={key:make_index(values) for key,values in raw_annotation.items()}
            for accession in ('SRR16194948','SRR22047451'):
                folder=args.out/accession;folder.mkdir();old=args.task/'rna_map_v1'/accession
                for item in mapping['libraries'][accession]['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
                bam=old/'star.Aligned.sortedByCoord.out.bam'
                command=[str(prefix/'bin/featureCounts'),'-T',str(os.environ['SLURM_CPUS_PER_TASK']),'-p','--countReadPairs',
                         '-s','0','-t','exon','-g','gene_id','-a',str(gtf),'-o',str(folder/'fragments.tsv'),
                         '-R','CORE','--Rpath',str(folder),str(bam)]
                with (folder/'featureCounts.stdout').open('w') as out,(folder/'featureCounts.stderr').open('w') as error:
                    subprocess.run(command,env=env,stdout=out,stderr=error,check=True)
                commands.append(command)
                assert count_table(folder/'fragments.tsv')==count_table(old/'fragments.tsv')
                summary=lambda path:[line.split('\t') for line in path.read_text().splitlines()[1:]]
                assert summary(folder/'fragments.tsv.summary')==summary(old/'fragments.tsv.summary')
                paths=list(folder.glob('*.featureCounts'));assert len(paths)==1
                core=defaultdict(list)
                with paths[0].open() as handle:
                    for line in handle:
                        fields=line.rstrip('\n').split('\t');assert len(fields)==4
                        core[fields[0]].append(fields[1:])
                primary=defaultdict(list)
                with (folder/'samtools.stderr').open('w') as error:
                    process=subprocess.Popen([str(prefix/'bin/samtools'),'view','-F','2304',str(bam)],env=env,
                                             stdout=subprocess.PIPE,stderr=error,text=True)
                    try:
                        for line in process.stdout:
                            fields=line.rstrip('\n').split('\t');primary[fields[0]].append(fields)
                        assert process.wait()==0
                    finally:
                        process.stdout.close()
                        if process.poll() is None:process.terminate();process.wait()
                eligibility=Counter();eligible=defaultdict(list)
                for name,records in primary.items():
                    assert len(records)==2
                    if any(int(row[1])&4 for row in records):eligibility['unmapped_mate']+=1;continue
                    nh=[next((int(tag[5:]) for tag in row[11:] if tag.startswith('NH:i:')),None) for row in records]
                    if nh!=[1,1]:eligibility['nonunique_mate']+=1;continue
                    if len(core.get(name,[]))!=1:eligibility['not_one_CORE_record']+=1;continue
                    status=core[name][0][0];eligible[status].append(name);eligibility['eligible_'+status]+=1
                selected=[]
                for status,n in [('Unassigned_NoFeatures',50),('Assigned',10),('Unassigned_Ambiguity',10)]:
                    selected.extend((name,status) for name in sorted(eligible[status])[:n])
                panel=[];counts=Counter();brute_checks=0
                for number,(name,status) in enumerate(selected):
                    genes={kind:set() for kind in ('gene','exon')};mates=[]
                    for row in primary[name]:
                        chrom=row[2];pieces=blocks(int(row[3]),row[5],len(row[9]));mate=dict(flag=int(row[1]),chromosome=chrom,
                            POS=int(row[3]),CIGAR=row[5],NH=1,sequence=row[9],quality=row[10],blocks=pieces)
                        for kind in ('gene','exon'):
                            key=(chrom,kind)
                            for start,end in pieces:
                                found=overlap(annotation[key],start,end) if key in annotation else set()
                                genes[kind].update(found)
                                if number<4:
                                    brute={g for a,b,g in raw_annotation.get(key,[]) if a<end and b>start}
                                    assert found==brute;brute_checks+=1
                        mates.append(mate)
                    label='exonic_support' if genes['exon'] else 'gene_body_only' if genes['gene'] else 'no_annotated_gene_body'
                    counts[status+'|'+label]+=1
                    reported=core[name][0][2].split(',') if status=='Assigned' else []
                    if reported and not set(reported)<=genes['exon']:counts['Assigned_reported_gene_without_exon_support']+=1
                    panel.append(dict(name=name,status=status,CORE=core[name],mates=mates,
                                      exon_genes=sorted(genes['exon']),gene_body_genes=sorted(genes['gene']),classification=label))
                panels[accession]=panel
                results[accession]=dict(eligibility=dict(eligibility),eligible_by_status={key:len(value) for key,value in eligible.items()},
                    selected_by_status=dict(Counter(row['status'] for row in panel)),classifications=dict(counts),
                    actual_annotation_brute_checks=brute_checks,original_gene_counts_and_summary_exact=True,
                    files={path.name:file_record(path) for path in sorted(folder.iterdir()) if path.is_file()})
            (args.out/'panels.json').write_text(json.dumps(panels,indent=2)+'\n')
            (args.out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
            (args.out/'report.json').write_text(json.dumps(dict(libraries=results,boundary_controls=controls,CIGAR_control=True,
                biological_null_specimens_verified=0,population_performance_estimated=False,count_default_changed=False),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
