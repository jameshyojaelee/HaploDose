"""Native DNA mapping of the two bounded same-library candidate read layers."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False);bwa=Path(os.environ['HAPLODOSE_BWA_BIN']).resolve()
    inputs=dict(script=Path(__file__),design=args.task/'DNA_MAP_DESIGN.md',bwa=bwa,
                wrapper=Path(__file__).resolve().parents[2]/'scripts/scone_bwa_readonly.sbatch',
                modules=Path(os.environ['HAPLODOSE_BWA_MODULE_RECORD']))
    for name in ('dna_index_v2','prepare_v1','prepare_audit_v1','tool_environment_v1'):
        inputs[name]=args.task/name/'run_manifest.json'
    run=start_run(args.out,inputs,dict(MAPQ_filter=None,UMI_deduplication=False),
                  dict(scope='Native DNA mapping only; no genomic CNA or clone truth'),mode='scone_dna_map')
    try:
        for name in ('dna_index_v2','prepare_v1','prepare_audit_v1','tool_environment_v1'):
            manifest=json.loads(inputs[name].read_text());assert manifest['status']=='complete'
            for item in manifest['output']['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
        index=json.loads((args.task/'dna_index_v2/report.json').read_text());runtime=index['native_runtime']
        assert file_record(bwa)['sha256']==runtime['executable']['sha256']
        assert os.environ['HAPLODOSE_BWA_LD_LIBRARY_PATH']==runtime['LD_LIBRARY_PATH']
        for item in runtime['libraries']:assert file_record(item['path'])['sha256']==item['sha256']
        tool=json.loads((args.task/'tool_environment_v1/report.json').read_text());prefix=Path(tool['environment'])
        for item in json.loads((args.task/'tool_environment_v1/tool_files.json').read_text()).values():
            assert file_record(item['path'])['sha256']==item['sha256']
        prepared=json.loads((args.task/'prepare_v1/report.json').read_text());outputs={};commands=[]
        for accession in ('SRR16194948','SRR22047451'):
            source=prepared['libraries'][accession];folder=args.out/accession;folder.mkdir()
            for item in source['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
            bwa_env={**os.environ,'LD_LIBRARY_PATH':runtime['LD_LIBRARY_PATH']}
            sam_env={**os.environ,'LD_LIBRARY_PATH':str(prefix/'lib')}
            rg='@RG\\tID:'+accession+'\\tSM:'+source['layout']['sample']
            command=[str(bwa),'mem','-t',str(os.environ['SLURM_CPUS_PER_TASK']),'-M','-R',rg,index['index_prefix'],
                     source['files']['DNA_R1.fastq.gz']['path'],source['files']['DNA_R2.fastq.gz']['path']]
            started=time.monotonic()
            with (folder/'alignments.sam').open('w') as out,(folder/'bwa.stderr').open('w') as error:
                subprocess.run(command,env=bwa_env,stdout=out,stderr=error,check=True)
            commands.append(dict(library=accession,label='bwa_mem',argv=command,seconds=time.monotonic()-started))
            for label,command in [
                ('view',[str(prefix/'bin/samtools'),'view','-b','-o',str(folder/'unsorted.bam'),str(folder/'alignments.sam')]),
                ('sort',[str(prefix/'bin/samtools'),'sort','-@','2','-m','512M','-o',str(folder/'DNA.bam'),str(folder/'unsorted.bam')]),
                ('flagstat',[str(prefix/'bin/samtools'),'flagstat',str(folder/'DNA.bam')])]:
                started=time.monotonic()
                with (folder/(label+'.stdout')).open('w') as out,(folder/(label+'.stderr')).open('w') as error:
                    subprocess.run(command,env=sam_env,stdout=out,stderr=error,check=True)
                commands.append(dict(library=accession,label=label,argv=command,seconds=time.monotonic()-started))
            (args.out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
            outputs[accession]=dict(input_candidate_DNA_pairs=source['counts']['DNA'],sample=source['layout']['sample'],
                files={p.name:file_record(p) for p in sorted(folder.iterdir()) if p.is_file()})
        (args.out/'report.json').write_text(json.dumps(dict(libraries=outputs,commands=commands,MAPQ_filter=None,
             PCR_UMI_deduplicated=False,DNA_copy_number_truth=False,biological_null_specimens_verified=0),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
