"""Independent metadata-only reconstruction, missing-value counterexample and gaps."""
import argparse
from collections import Counter,defaultdict
import csv
from decimal import Decimal,InvalidOperation
import importlib.util
import json
from pathlib import Path
import re

from haplodose.provenance import start_run,finish_run,file_record


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def fixtures(current,baseline):
    def record(name,cell,layer,**changes):
        characteristics={'cell id':cell,'cell line':'test_line','cell type':'test_type','cell count':'1','treatment':'untreated','time point':'0h'}
        characteristics.update(changes)
        return dict(accession=name,fields={'molecule_ch1':['total RNA' if layer=='RNA' else 'genomic DNA']},characteristics=characteristics)
    unknown=[record('RNA','A','RNA',treatment='NA'),record('DNA','A','DNA',treatment='NA')]
    assert baseline.pairs(unknown)[0][0]['status']=='paired_metadata_only'
    assert current.pairs(unknown)[0][0]['status']=='discordant_or_missing_metadata'
    rows=[record('R1','A','RNA'),record('D1','A','DNA'),record('R2','B','RNA'),record('D2','B','DNA',treatment='treated'),
          record('R3','C','RNA',**{'cell count':'2'}),record('D3','C','DNA',**{'cell count':'2'}),
          record('R4','D','RNA'),record('R5','E','RNA'),record('R6','E','RNA'),record('D5','E','DNA')]
    expected={'A':'paired_metadata_only','B':'discordant_or_missing_metadata','C':'not_certified_single_cell',
              'D':'unresolved_library_cardinality','E':'unresolved_library_cardinality'}
    for ordered in (rows,list(reversed(rows))):
        pairs,unkeyed=current.pairs(ordered);assert not unkeyed
        assert {x['cell_id']:x['status'] for x in pairs}==expected
    noid=record('R7','F','RNA');noid['characteristics'].pop('cell id')
    assert current.pairs([noid])==([],['R7'])
    return dict(status_fixtures=11,baseline_NA_counterexample_reproduced=True,corrected_NA_rejected=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);root=a.task.parents[1];source=a.task/'metadata_v2'
    manifest=json.loads((source/'run_manifest.json').read_text());assert manifest['status']=='complete'
    for group in ('input','output'):
        for record in manifest[group]['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
    script=Path(__file__).with_name('dntr_correspondence_metadata.py');baseline=a.task/'source_v1/benchmarks/hardening/dntr_correspondence_metadata.py'
    assert file_record(script)['sha256']==manifest['input']['files']['script']['sha256']
    inputs=dict(script=Path(__file__),candidate=script,baseline=baseline,design=a.task/'DESIGN.md',
        missing_value_amendment=a.task/'MISSING_VALUE_AMENDMENT.md',metadata=source/'run_manifest.json')
    for name in ('GEO_metadata_only.soft','GEO_samples.json','ENA_runs.tsv','run_correspondence.json','cell_correspondence.json','file_presence.json','report.json'):
        inputs[name]=source/name
    inputs['frozen_runs']=root/'benchmarks/heldout2/dntrseq_gse144296_runs.tsv'
    inputs['download_manifest']=root/'data/external/dntrseq_gse144296/manifest.tsv'
    rec=start_run(a.out,inputs,{},dict(scope='Administrative library correspondence,not molecular match or independent biological truth'),mode='dntr_correspondence_audit')
    try:
        current=load(script,'corrected_dntr_metadata');old=load(baseline,'initial_dntr_metadata');check=fixtures(current,old)
        # Parse the whitelisted source without calling the candidate parser.
        rebuilt=[]
        for block in inputs['GEO_metadata_only.soft'].read_text().split('^SAMPLE = ')[1:]:
            lines=block.splitlines();row=dict(accession=lines[0],fields={},characteristics={})
            for line in lines[1:]:
                key,value=line.split(' = ',1);key=key[len('!Sample_'):]
                if key=='characteristics_ch1':
                    name,content=value.split(': ',1);name=name.lower().strip()
                    assert name not in row['characteristics'];row['characteristics'][name]=content
                else:row['fields'].setdefault(key,[]).append(value)
            rebuilt.append(row)
        records=json.loads(inputs['GEO_samples.json'].read_text());assert rebuilt==records
        assert len({r['accession'] for r in records})==len(records)
        by_sample={r['accession']:r for r in records};by_exp=defaultdict(set)
        assemblies=Counter();processing=Counter()
        for row in records:
            for text in row['fields'].get('data_processing',[]):
                processing[text]+=1
                for assembly in re.findall(r'Genome_build:\s*(\S+)',text):assemblies[assembly]+=1
            for relation in row['fields'].get('relation',[]):
                for accession in re.findall(r'SRX[0-9]+',relation):by_exp[accession].add(row['accession'])
        ena=list(csv.DictReader(inputs['ENA_runs.tsv'].open(),delimiter='\t'));ena_by_run={r['run_accession']:r for r in ena}
        frozen=list(csv.DictReader(inputs['frozen_runs'].open(),delimiter='\t'));old_by_run={r['run']:r for r in frozen}
        runs=json.loads(inputs['run_correspondence.json'].read_text());assert len(runs)==len(old_by_run)==4869
        assert set(ena_by_run)==set(old_by_run)
        layer_conflicts=[];line_conflicts=[];run_count=Counter()
        for row in runs:
            new,previous=ena_by_run[row['run']],old_by_run[row['run']]
            accessions=by_exp[new['experiment_accession']];assert set(row['samples'])==accessions and len(accessions)==1
            sample=by_sample[next(iter(accessions))];molecule=sample['fields']['molecule_ch1']
            layer='DNA' if molecule==['genomic DNA'] else 'RNA' if molecule==['total RNA'] else 'unresolved'
            assert row['roles']==[layer] and row['explicit_cell_ids']==[sample['characteristics']['cell id']]
            if layer.lower()!=previous['group'].split('_')[-1]:layer_conflicts.append(row['run'])
            if sample['characteristics']['cell line']!=previous['group'].split('_')[-2]:line_conflicts.append(row['run'])
            assert row['read_metadata_unchanged']==all(new[k]==previous[k] for k in ('read_count','base_count','fastq_md5','fastq_bytes'))
            run_count[sample['accession']]+=1
        by_cell=defaultdict(list)
        for row in records:by_cell[row['characteristics']['cell id']].append(row)
        cell_rows=json.loads(inputs['cell_correspondence.json'].read_text());assert set(by_cell)=={x['cell_id'] for x in cell_rows}
        statuses=Counter();treatments=Counter()
        missing={'','na','n/a','unknown','not available','not provided','not applicable'}
        for row in cell_rows:
            samples=by_cell[row['cell_id']];dna=[s for s in samples if s['fields']['molecule_ch1']==['genomic DNA']]
            rna=[s for s in samples if s['fields']['molecule_ch1']==['total RNA']]
            assert set(row['dna_samples'])=={s['accession'] for s in dna};assert set(row['rna_samples'])=={s['accession'] for s in rna}
            if len(dna)!=1 or len(rna)!=1 or len(samples)!=2:status='unresolved_library_cardinality'
            else:
                values=[s['characteristics'] for s in samples];known_equal=True
                for name in ('cell line','cell type','cell count','treatment','time point'):
                    pair=[v.get(name) for v in values]
                    if pair[0]!=pair[1] or any(v is None or str(v).strip().lower() in missing for v in pair):known_equal=False
                try:single=all(Decimal(v.get('cell count',''))==1 for v in values)
                except InvalidOperation:single=False
                status='not_certified_single_cell' if not single else 'paired_metadata_only' if known_equal else 'discordant_or_missing_metadata'
            assert row['status']==status,(row['cell_id'],row['status'],status);statuses[status]+=1
            annotation=samples[0]['characteristics'];treatments[annotation['cell line'],annotation['treatment'],annotation['time point'],status]+=1
        archived=list(csv.DictReader(inputs['download_manifest'].open(),delimiter='\t'));by_file={r['file']:r for r in archived}
        files=json.loads(inputs['file_presence.json'].read_text());gaps=[]
        for row in files:
            new=ena_by_run[row['run']];basename=Path(row['file']).name;names=new['fastq_ftp'].split(';')
            index=[Path(s).name for s in names].index(basename);expected_md5=new['fastq_md5'].split(';')[index];expected_size=new['fastq_bytes'].split(';')[index]
            old=by_file.get(row['file']);matches=bool(old and old['md5']==expected_md5 and old['bytes']==expected_size)
            assert row['archived_metadata_matches']==matches
            path=root/'data/external/dntrseq_gse144296'/row['file'];assert row['exists']==path.is_file()
            assert row['size_matches']==(path.is_file() and path.stat().st_size==int(expected_size))
            if not row['size_matches']:
                gaps.append(dict(**row,url='https://'+names[index],expected_md5=expected_md5,expected_bytes=int(expected_size)))
        report=json.loads(inputs['report.json'].read_text());assert dict(statuses)==report['cell_status']
        assert len(files)==report['expected_files'] and len(gaps)==report['missing_or_wrong_size']
        assert not report['cohort_admitted'] and report['certified_single_CNA_null_specimens']==0
        (a.out/'missing_files.json').write_text(json.dumps(gaps,indent=2)+'\n')
        (a.out/'treatment_status.json').write_text(json.dumps([dict(cell_line=k[0],treatment=k[1],time_point=k[2],status=k[3],cell_ids=n) for k,n in sorted(treatments.items())],indent=2)+'\n')
        (a.out/'report.json').write_text(json.dumps(dict(all_checks_passed=True,fixtures=check,geo_records_reconstructed=len(rebuilt),run_joins_checked=len(runs),
            cell_ids_checked=len(cell_rows),cell_status=dict(statuses),file_metadata_checks=len(files),missing_files=len(gaps),
            missing_bytes=sum(r['expected_bytes'] for r in gaps),layer_conflicts=layer_conflicts,line_conflicts=line_conflicts,
            largest_technical_runs_per_GEO_library=max(run_count.values()),biological_null_specimens_certified=0,
            declared_assembly_records=dict(assemblies),data_processing_records=dict(processing),
            raw_reads_opened=False,measurements_opened=False),indent=2)+'\n')
        for record in rec['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(a.out,rec,backend='CPU independent metadata parsing,constraint joins and file-stat checks')
    except BaseException as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
