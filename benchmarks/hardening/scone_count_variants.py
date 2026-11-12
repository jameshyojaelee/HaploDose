"""Fixed matched RNA fragment-count definitions; no scientific default promotion."""
import argparse
from collections import Counter, defaultdict
import csv
import fcntl
import gzip
import json
import os
from pathlib import Path
import subprocess

from haplodose.provenance import start_run, finish_run, file_record


def table(path):
    with path.open() as handle:
        rows=list(csv.reader((line for line in handle if not line.startswith('#')),delimiter='\t'))
    counts={row[0]:int(row[-1]) for row in rows[1:]};assert len(counts)==len(rows)-1 and min(counts.values())>=0
    return counts


def core(path):
    records=defaultdict(list)
    with path.open() as handle:
        for line in handle:
            fields=line.rstrip('\n').split('\t');assert len(fields)==4
            records[fields[0]].append(fields[1:])
    return records


def assigned(records):
    names={};counts=Counter()
    for name,rows in records.items():
        if any(row[0]=='Assigned' for row in rows):
            assert len(rows)==1 and rows[0][0]=='Assigned' and rows[0][1]=='1'
            target=rows[0][2];assert target not in ('','NA') and ',' not in target
            names[name]=target;counts[target]+=1
    return names,counts


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    inputs=dict(script=Path(__file__),design=args.task/'COUNT_VARIANTS_DESIGN.md')
    for name in ('nofeatures_v1','nofeatures_audit_v1','rna_map_v1','reference_v3','tool_environment_v1'):
        inputs[name]=args.task/name/'run_manifest.json'
    run=start_run(args.out,inputs,dict(definitions=['exon','gene_body','exon_priority_rescue']),
                  dict(scope='Matched measurement definitions; no calibration or DNA accuracy claim'),mode='scone_count_variants')
    try:
        for name,path in inputs.items():
            if name in ('script','design'):continue
            manifest=json.loads(path.read_text());assert manifest['status']=='complete'
            for item in manifest['output']['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
        tools=json.loads((args.task/'tool_environment_v1/report.json').read_text());prefix=Path(tools['environment'])
        inventory=json.loads((args.task/'tool_environment_v1/tool_files.json').read_text())
        for item in inventory.values():assert file_record(item['path'])['sha256']==item['sha256']
        env={**os.environ,'LD_LIBRARY_PATH':str(prefix/'lib')}
        reference=json.loads((args.task/'reference_v3/report.json').read_text())
        gtf=reference['reference']['Homo_sapiens.GRCh38.112.gtf.gz']
        old_report=json.loads((args.task/'nofeatures_v1/report.json').read_text())
        mapped=json.loads((args.task/'rna_map_v1/report.json').read_text());results={};commands=[]
        with (Path(gtf['path']).parent.parent/'.setup_reference.lock').open('rb') as lock:
            fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB);assert file_record(gtf['path'])['sha256']==gtf['sha256']
            for accession in ('SRR16194948','SRR22047451'):
                for report in (old_report,mapped):
                    for item in report['libraries'][accession]['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
                folder=args.out/accession;folder.mkdir();old=args.task/'nofeatures_v1'/accession
                bam=args.task/'rna_map_v1'/accession/'star.Aligned.sortedByCoord.out.bam'
                command=[str(prefix/'bin/featureCounts'),'-T',str(os.environ['SLURM_CPUS_PER_TASK']),'-p','--countReadPairs',
                     '-s','0','-t','gene','-g','gene_id','-a',gtf['path'],'-o',str(folder/'gene_body.tsv'),
                     '-R','CORE','--Rpath',str(folder),str(bam)]
                with (folder/'featureCounts.stdout').open('w') as out,(folder/'featureCounts.stderr').open('w') as error:
                    subprocess.run(command,env=env,stdout=out,stderr=error,check=True)
                commands.append(command)
                paths=list(folder.glob('*.featureCounts'));assert len(paths)==1
                body_core=core(paths[0]);exon_core=core(next(old.glob('*.featureCounts')))
                assert set(body_core)==set(exon_core)
                body_names,body_tally=assigned(body_core);exon_names,exon_tally=assigned(exon_core)
                body=table(folder/'gene_body.tsv');exon=table(old/'fragments.tsv')
                assert {gene:value for gene,value in body.items() if value}==dict(body_tally)
                assert {gene:value for gene,value in exon.items() if value}==dict(exon_tally)
                assert sum(exon.values())==mapped['libraries'][accession]['exonic_assigned_fragments']
                rescue={name:gene for name,gene in body_names.items()
                    if len(exon_core[name])==1 and exon_core[name][0][0]=='Unassigned_NoFeatures'}
                assert not set(rescue)&set(exon_names)
                tally=Counter(rescue.values());genes=sorted(set(exon)|set(body));priority={g:exon.get(g,0)+tally[g] for g in genes}
                assert all(priority[g]>=exon.get(g,0) for g in genes)
                assert sum(priority.values())==sum(exon.values())+len(rescue)
                with (folder/'count_definitions.tsv').open('w') as handle:
                    writer=csv.writer(handle,delimiter='\t');writer.writerow(['gene_id','exon','gene_body','exon_priority_rescue'])
                    writer.writerows((g,exon.get(g,0),body.get(g,0),priority[g]) for g in genes)
                with gzip.open(folder/'rescued_templates.tsv.gz','wt') as handle:
                    writer=csv.writer(handle,delimiter='\t');writer.writerow(['QNAME_candidate_UMI','gene_id'])
                    writer.writerows(sorted(rescue.items()))
                summary={line.split('\t')[0]:int(line.split('\t')[1]) for line in (folder/'gene_body.tsv.summary').read_text().splitlines()[1:]}
                assert summary['Assigned']==sum(body.values())
                results[accession]=dict(gene_universe=len(genes),exonic_assigned=sum(exon.values()),gene_body_assigned=sum(body.values()),
                    exon_priority_assigned=sum(priority.values()),rescued_templates=len(rescue),
                    exon_assigned_body_not_assigned=len(set(exon_names)-set(body_names)),
                    exon_body_disagreeing_target=sum(body_names.get(name,gene)!=gene for name,gene in exon_names.items() if name in body_names),
                    gene_body_summary=summary,CORE_native_gene_tallies_exact=True,rescue_conservation_exact=True,
                    files={p.name:file_record(p) for p in sorted(folder.iterdir()) if p.is_file()})
            (args.out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
            (args.out/'report.json').write_text(json.dumps(dict(libraries=results,count_unit='PCR-affected RNA read fragments',
                biological_null_specimens_verified=0,raw_UMI_molecules=False,count_default_changed=False,
                scientific_accuracy_evaluated=False),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
