"""Independent one-based, observed-position audit of the entire frozen read panel."""
import argparse
from bisect import bisect_left
from collections import Counter, defaultdict
import gzip
import json
import os
from pathlib import Path
import subprocess

from haplodose.provenance import start_run, finish_run, file_record


def observed_positions(pos,cigar,query_length):
    current=pos;observed=[];query=0;digits=''
    for char in cigar:
        if char.isdigit():digits+=char;continue
        assert digits;length=int(digits);digits='';assert length>0
        if char in ('M','=','X'):
            observed.extend(range(current,current+length));current+=length;query+=length
        elif char in ('N','D'):current+=length
        elif char in ('S','I'):query+=length
        else:assert char in ('H','P')
    assert not digits and query==query_length
    return observed


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--prepared-read-audit',action='store_true')
    args=parser.parse_args();args.out.mkdir(exist_ok=False);parent=args.task/'nofeatures_v1'
    design='NOFEATURES_RAW_AUDIT_DESIGN.md' if args.prepared_read_audit else 'NOFEATURES_AUDIT_DESIGN.md'
    inputs=dict(script=Path(__file__),design=args.task/design,original=parent/'run_manifest.json')
    if args.prepared_read_audit:inputs['prepared']=args.task/'prepare_v1/run_manifest.json'
    run=start_run(args.out,inputs,dict(prepared_read_audit=args.prepared_read_audit),
                  dict(scope='Exact computation replication; no new biological hypothesis or selected panel'),mode='scone_nofeatures_audit')
    try:
        manifest=json.loads((parent/'run_manifest.json').read_text());assert manifest['status']=='complete'
        for section in ('input','output'):
            for item in manifest[section]['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
        report=json.loads((parent/'report.json').read_text());panels=json.loads((parent/'panels.json').read_text())
        tools=json.loads((args.task/'tool_environment_v1/report.json').read_text());prefix=Path(tools['environment'])
        assert file_record(prefix/'bin/samtools')['sha256']==tools['versions']['samtools']['executable']['sha256']
        mapping=json.loads((args.task/'rna_map_v1/report.json').read_text())
        prepared=None
        if args.prepared_read_audit:
            prepared_manifest=json.loads(inputs['prepared'].read_text());assert prepared_manifest['status']=='complete'
            for item in prepared_manifest['output']['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
            prepared=json.loads((args.task/'prepare_v1/report.json').read_text())
        env={**os.environ,'LD_LIBRARY_PATH':str(prefix/'lib')};indexed=defaultdict(list);reference_sets={};checks={}
        for accession,entries in panels.items():
            original=report['libraries'][accession]
            for item in original['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
            core=defaultdict(list)
            path=next(parent.joinpath(accession).glob('*.featureCounts'))
            for line in path.read_text().splitlines():
                name,*fields=line.split('\t');assert len(fields)==3;core[name].append(fields)
            bam=args.task/'rna_map_v1'/accession/'star.Aligned.sortedByCoord.out.bam'
            for record in mapping['libraries'][accession]['files'].values():
                assert file_record(record['path'])['sha256']==record['sha256']
            primary=defaultdict(list)
            with (args.out/(accession+'.samtools.stderr')).open('w') as error:
                process=subprocess.Popen([str(prefix/'bin/samtools'),'view','-F','2304',str(bam)],env=env,
                    stdout=subprocess.PIPE,stderr=error,text=True)
                try:
                    for line in process.stdout:
                        fields=line.rstrip('\n').split('\t');primary[fields[0]].append(fields)
                    assert process.wait()==0
                finally:
                    process.stdout.close()
                    if process.poll() is None:process.terminate();process.wait()
            eligible=defaultdict(list);eligibility=Counter()
            for name,records in primary.items():
                assert len(records)==2
                if any(int(row[1])&4 for row in records):eligibility['unmapped_mate']+=1;continue
                if any('NH:i:1' not in row[11:] for row in records):eligibility['nonunique_mate']+=1;continue
                if len(core.get(name,[]))!=1:eligibility['not_one_CORE_record']+=1;continue
                status=core[name][0][0];eligible[status].append(name);eligibility['eligible_'+status]+=1
            expected=[(name,status) for status,n in [('Unassigned_NoFeatures',50),('Assigned',10),('Unassigned_Ambiguity',10)]
                      for name in sorted(eligible[status])[:n]]
            assert expected==[(item['name'],item['status']) for item in entries]
            assert dict(eligibility)==original['eligibility']
            original_reads={}
            if prepared is not None:
                selected_names={item['name'] for item in entries}
                for mate in (1,2):
                    record=prepared['libraries'][accession]['files'][f'RNA_R{mate}.fastq.gz']
                    assert file_record(record['path'])['sha256']==record['sha256']
                    with gzip.open(record['path'],'rt',encoding='ascii') as handle:
                        while True:
                            header=handle.readline()
                            if not header:break
                            sequence=handle.readline().rstrip('\r\n');plus=handle.readline();quality=handle.readline().rstrip('\r\n')
                            assert header.startswith('@') and plus.startswith('+') and len(sequence)==len(quality)
                            name=header[1:].split()[0]
                            if name in selected_names:
                                key=(name,mate);assert key not in original_reads;original_reads[key]=(sequence,quality)
                assert len(original_reads)==2*len(entries)
            for item in entries:
                name=item['name'];key=(accession,name);reference_sets[key]={'gene':set(),'exon':set()}
                assert item['CORE']==core[name]
                records={1 if int(row[1])&64 else 2:row for row in primary[name]}
                for mate in item['mates']:
                    row=records[1 if mate['flag']&64 else 2]
                    assert [mate['chromosome'],str(mate['POS']),mate['CIGAR'],mate['sequence'],mate['quality']]==[row[i] for i in (2,3,5,9,10)]
                    if prepared is not None:
                        sequence,quality=row[9:11]
                        if int(row[1])&16:
                            sequence=sequence.translate(str.maketrans('ACGTNacgtn','TGCANtgcan'))[::-1];quality=quality[::-1]
                        raw_mate=2 if int(row[1])&64 else 1
                        assert (sequence,quality)==original_reads[name,raw_mate],f'Raw prepared mate mismatch: {name}/{raw_mate}'
                    positions=observed_positions(int(row[3]),row[5],len(row[9]))
                    assert positions==[p+1 for a,b in mate['blocks'] for p in range(a,b)]
                    if positions:indexed[row[2]].append((key,positions))
            checks[accession]=dict(selected=len(entries),selected_identities_exact=True,primary_read_correspondence_exact=True)
            if prepared is not None:checks[accession]['prepared_FASTQ_mate_sequences_qualities_exact']=len(original_reads)
        ref=json.loads((args.task/'reference_v3/report.json').read_text())['reference']['Homo_sapiens.GRCh38.112.gtf.gz']
        assert file_record(ref['path'])['sha256']==ref['sha256']
        with gzip.open(ref['path'],'rt') as handle:
            for line in handle:
                if line.startswith('#'):continue
                fields=line.rstrip('\n').split('\t');chrom,kind=fields[0],fields[2]
                if kind not in ('gene','exon') or chrom not in indexed:continue
                start,end=int(fields[3]),int(fields[4]);attributes={}
                for part in fields[8].split(';'):
                    parts=part.strip().split(' ',1)
                    if len(parts)==2:attributes[parts[0]]=parts[1].strip('"')
                gene=attributes['gene_id']
                for key,positions in indexed[chrom]:
                    if end<positions[0] or start>positions[-1]:continue
                    offset=bisect_left(positions,start)
                    if offset<len(positions) and positions[offset]<=end:reference_sets[key][kind].add(gene)
        totals={}
        for accession,entries in panels.items():
            classified=Counter();ambiguity=Counter();assigned_checked=0
            for item in entries:
                sets=reference_sets[accession,item['name']]
                assert sets['gene']==set(item['gene_body_genes']) and sets['exon']==set(item['exon_genes'])
                label='exonic_support' if sets['exon'] else 'gene_body_only' if sets['gene'] else 'no_annotated_gene_body'
                assert label==item['classification'];classified[item['status']+'|'+label]+=1
                if item['status']=='Assigned':
                    assert set(item['CORE'][0][2].split(','))<=sets['exon'];assigned_checked+=1
                elif item['status']=='Unassigned_Ambiguity':
                    ambiguity['multiple_exonic_genes' if len(sets['exon'])>=2 else 'unresolved_single_or_no_gene']+=1
            assert dict(classified)==report['libraries'][accession]['classifications']
            totals[accession]={**checks[accession], 'classifications':dict(classified),
                'assigned_target_genes_exact':assigned_checked,'ambiguity_gene_support':dict(ambiguity),
                'all_observed_positions_and_gene_sets_exact':True}
        (args.out/'report.json').write_text(json.dumps(dict(libraries=totals,all_selected_templates=sum(len(v) for v in panels.values()),
            one_based_direct_annotation_reference=True,biological_null_specimens_verified=0,
            population_validation=False,count_default_changed=False),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
