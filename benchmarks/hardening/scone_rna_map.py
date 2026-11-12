"""Declared RNA mapping and read-fragment quantification for the two bounded libraries."""
import argparse
import csv
import fcntl
import json
import os
from pathlib import Path
import subprocess
import time

from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    inputs=dict(script=Path(__file__),design=args.task/'RNA_MAP_DESIGN.md')
    for name in ('reference_v3','prepare_audit_v1','prepare_v1','tool_environment_v1'):
        inputs[name]=args.task/name/'run_manifest.json'
    run=start_run(args.out,inputs,dict(count_unit='PCR-affected exonic RNA read fragments',UMI_deduplication=False),
                  dict(scope='Two development read subsets; no CNA inference or biological validation'),mode='scone_rna_map')
    try:
        for name,path in inputs.items():
            if name in ('script','design'):continue
            manifest=json.loads(path.read_text());assert manifest['status']=='complete'
            for record in manifest['output']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        reference=json.loads((args.task/'reference_v3/report.json').read_text())
        prepared=json.loads((args.task/'prepare_v1/report.json').read_text())
        tool_report=json.loads((args.task/'tool_environment_v1/report.json').read_text())
        tools=Path(tool_report['environment']);env={**os.environ,'LD_LIBRARY_PATH':str(tools/'lib')}
        inventory=json.loads((args.task/'tool_environment_v1/tool_files.json').read_text())
        for record in inventory.values():assert file_record(record['path'])['sha256']==record['sha256']
        refdir=Path(reference['reference']['Homo_sapiens.GRCh38.dna.primary_assembly.fa']['path']).parent
        outputs={};commands=[]
        with (refdir.parent/'.setup_reference.lock').open('rb') as lock:
            fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
            reference_files=[*reference['reference'].values(),*reference['index_artifacts']]
            for record in reference_files:assert file_record(record['path'])['sha256']==record['sha256']
            stamps={r['path']:(Path(r['path']).stat().st_size,Path(r['path']).stat().st_mtime_ns) for r in reference_files}
            for accession in ('SRR16194948','SRR22047451'):
                entry=prepared['libraries'][accession]
                for item in entry['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
                folder=args.out/accession;folder.mkdir()
                first=entry['files']['RNA_R2.fastq.gz']['path'];second=entry['files']['RNA_R1.fastq.gz']['path']
                tasks=[('STAR',[str(tools/'bin/STAR-plain'),'--runThreadN',str(os.environ['SLURM_CPUS_PER_TASK']),
                   '--genomeDir',str(refdir/'star_index'),'--genomeLoad','NoSharedMemory','--readFilesIn',first,second,
                   '--readFilesCommand','/usr/bin/gzip','-cd','--outFilterType','BySJout','--outFilterMultimapNmax','20',
                   '--alignSJoverhangMin','8','--alignSJDBoverhangMin','1','--outFilterMismatchNmax','999',
                   '--outFilterMismatchNoverReadLmax','0.04','--alignIntronMin','20','--alignIntronMax','1000000',
                   '--alignMatesGapMax','1000000','--outSAMtype','BAM','SortedByCoordinate','--limitBAMsortRAM','8000000000',
                   '--outSAMattributes','NH','HI','AS','nM','--outSAMunmapped','Within','KeepPairs',
                   '--quantMode','GeneCounts','--outFileNamePrefix',str(folder/'star.')]),
                  ('flagstat',[str(tools/'bin/samtools'),'flagstat',str(folder/'star.Aligned.sortedByCoord.out.bam')]),
                  ('featureCounts',[str(tools/'bin/featureCounts'),'-T',str(os.environ['SLURM_CPUS_PER_TASK']),
                   '-p','--countReadPairs','-s','0','-t','exon','-g','gene_id','-a',
                   str(refdir/'Homo_sapiens.GRCh38.112.gtf.gz'),'-o',str(folder/'fragments.tsv'),
                   str(folder/'star.Aligned.sortedByCoord.out.bam')])]
                for label,command in tasks:
                    started=time.monotonic()
                    with (folder/(label+'.stdout')).open('w') as out,(folder/(label+'.stderr')).open('w') as error:
                        subprocess.run(command,env=env,stdout=out,stderr=error,check=True)
                    commands.append(dict(library=accession,label=label,argv=command,seconds=time.monotonic()-started))
                    (args.out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
                with (folder/'fragments.tsv').open() as handle:
                    rows=list(csv.DictReader((line for line in handle if not line.startswith('#')),delimiter='\t'))
                assert rows and len(rows[0])==7;column=list(rows[0])[-1]
                values=[int(row[column]) for row in rows];assert all(value>=0 for value in values)
                with (folder/'fragments.tsv.summary').open() as handle:
                    summary={row[0]:int(row[1]) for row in list(csv.reader(handle,delimiter='\t'))[1:]}
                assert sum(values)==summary['Assigned']
                outputs[accession]=dict(input_RNA_pairs=entry['counts']['RNA'],genes=len(values),
                   exonic_assigned_fragments=sum(values),featureCounts_summary=summary,
                   files={str(path.relative_to(args.out)):file_record(path) for path in sorted(folder.rglob('*')) if path.is_file()})
            assert all((Path(path).stat().st_size,Path(path).stat().st_mtime_ns)==stamp for path,stamp in stamps.items())
        (args.out/'report.json').write_text(json.dumps(dict(libraries=outputs,commands=commands,
             tools=tool_report['versions'],count_unit='PCR-affected exonic RNA read fragments',raw_UMI_molecules=False,
             calibration_validated=False,biological_null_specimens_verified=0),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
