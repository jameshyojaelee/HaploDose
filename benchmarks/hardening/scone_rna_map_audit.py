"""Independent primary-template and count-unit audit of the bounded RNA mappings."""
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
    args=parser.parse_args();args.out.mkdir(exist_ok=False);mapped=args.task/'rna_map_v1'
    run=start_run(args.out,dict(script=Path(__file__),design=args.task/'RNA_MAP_DESIGN.md',
                            mapping=mapped/'run_manifest.json',preparation=args.task/'prepare_v1/run_manifest.json'),{},
                  dict(scope='Primary pair conservation and output units, not independent gene truth'),mode='scone_rna_map_audit')
    try:
        for directory in (mapped,args.task/'prepare_v1'):
            manifest=json.loads((directory/'run_manifest.json').read_text());assert manifest['status']=='complete'
            for item in manifest['output']['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
        source=json.loads((args.task/'prepare_v1/report.json').read_text());result=json.loads((mapped/'report.json').read_text())
        tool=json.loads((args.task/'tool_environment_v1/report.json').read_text());prefix=Path(tool['environment'])
        env={**os.environ,'LD_LIBRARY_PATH':str(prefix/'lib')};audited={}
        for accession in ('SRR16194948','SRR22047451'):
            entry=source['libraries'][accession];expected=set()
            for item in entry['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
            with gzip.open(entry['files']['RNA_R2.fastq.gz']['path'],'rt') as handle:
                for index,line in enumerate(handle):
                    if index%4==0:
                        name=line.rstrip()[1:];assert name not in expected;expected.add(name)
            assert len(expected)==entry['counts']['RNA']
            for item in result['libraries'][accession]['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
            bam=mapped/accession/'star.Aligned.sortedByCoord.out.bam';seen=Counter();states=Counter()
            command=[str(prefix/'bin/samtools'),'view','-F','2304',str(bam)]
            with (args.out/(accession+'.samtools.stderr')).open('w') as error:
                process=subprocess.Popen(command,env=env,stdout=subprocess.PIPE,stderr=error,text=True)
                try:
                    for line in process.stdout:
                        fields=line.rstrip('\n').split('\t');assert len(fields)>=11
                        name,flag=fields[0],int(fields[1]);assert name in expected and flag&1
                        assert bool(flag&64)!=bool(flag&128)
                        mate=1 if flag&64 else 2;seen[name,mate]+=1
                        assert seen[name,mate]==1,'Repeated primary mate'
                        states['unmapped' if flag&4 else 'mapped']+=1
                    assert process.wait()==0
                finally:
                    process.stdout.close()
                    if process.poll() is None:process.terminate();process.wait()
            assert set(seen)=={(name,mate) for name in expected for mate in (1,2)}
            assert sum(seen.values())==2*len(expected)
            lines=[line.rstrip('\n').split('\t') for line in (mapped/accession/'fragments.tsv').read_text().splitlines() if not line.startswith('#')]
            assert lines[0][:6]==['Geneid','Chr','Start','End','Strand','Length']
            ids=[row[0] for row in lines[1:]];assert len(ids)==len(set(ids))
            values=[int(row[-1]) for row in lines[1:]];assert all(value>=0 for value in values)
            summary_lines=(mapped/accession/'fragments.tsv.summary').read_text().splitlines()
            summary={line.split('\t')[0]:int(line.split('\t')[1]) for line in summary_lines[1:]}
            assert sum(values)==summary['Assigned']==result['libraries'][accession]['exonic_assigned_fragments']
            audited[accession]=dict(primary_templates=len(expected),primary_records=sum(seen.values()),primary_states=dict(states),
                all_ID_UMI_mate_pairs_exact=True,exonic_assigned_fragments=sum(values),featureCounts_summary=summary,
                independently_validated_gene_assignments=False)
        (args.out/'report.json').write_text(json.dumps(dict(libraries=audited,count_unit='PCR-affected exonic RNA read fragments',
              raw_UMI_molecules=False,biological_null_specimens_verified=0),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
