"""Independently audit raw SNP roles, observation conservation and public outputs."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def inventory(path,required_role=None):
    rows=Counter();totals=[0,0,0];sites={};examples={}
    with path.open(newline='') as stream:
        reader=csv.DictReader(stream,delimiter='\t');fields=reader.fieldnames
        for row in reader:
            chrom=row['CHROM'] if row['CHROM'].startswith('chr') else 'chr'+row['CHROM']
            pos=int(row['POS']);key=chrom+':'+str(pos)
            payload='\0'.join(['haplodose-allele-role-v1','20261108','hg38',chrom,str(pos)])
            number=int(hashlib.blake2b(payload.encode()).hexdigest()[:16],16)
            role='discovery' if number<9223372036854775808 else 'reserved'
            if required_role:
                assert role==required_role
            sites[key]=role;examples.setdefault(role,row)
            rows[tuple(row[name] for name in fields)]+=1
            totals[0]+=1;totals[1]+=int(row['AD']);totals[2]+=int(row['DP'])
    return fields,rows,totals,sites,examples


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    source=args.task/'prepare_v1'
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md')
    manifests={}
    for label in ('prepare','binomial','beta'):
        directory=source if label=='prepare' else args.task/(label+'_v1')
        inputs[label]=directory/'run_manifest.json'
        manifest=json.loads(inputs[label].read_text());manifests[label]=manifest
        assert manifest['status']=='complete'
        for group in ('input','output'):
            for record in manifest[group]['files'].values():
                assert file_record(record['path'])['sha256']==record['sha256']
    helper=Path(manifests['prepare']['input']['files']['helper']['path']);inputs['helper']=helper
    spec=importlib.util.spec_from_file_location('frozen_roles',helper);roles=importlib.util.module_from_spec(spec);spec.loader.exec_module(roles)
    run=start_run(args.out,inputs,{},dict(scope='Role/IO/discovery identity only; no biological calibration'),
                  mode='allele_raw_role_audit')
    try:
        cases=json.loads((source/'report.json').read_text())['cases'];conservation=[];original=None
        for case in cases:
            fields,raw,totals,site_map,examples=inventory(source/(case+'.tsv'))
            dr=inventory(source/(case+'_roles/discovery.tsv'),'discovery')
            rr=inventory(source/(case+'_roles/reserved.tsv'),'reserved')
            assert fields==dr[0]==rr[0] and raw==dr[1]+rr[1]
            assert totals==[a+b for a,b in zip(dr[2],rr[2])]
            assert not set(dr[3])&set(rr[3]) and site_map==dr[3]|rr[3]
            declared=json.loads((source/(case+'_roles/roles.json')).read_text())
            assert declared['sites']==site_map
            if original is None:
                original=site_map;baseline_fields=fields;baseline_examples=examples
            assert all(site_map[key]==value for key,value in original.items())
            conservation.append(dict(case=case,input_rows=totals[0],raw_AD=totals[1],raw_DP=totals[2],
                                     discovery_rows=dr[2][0],reserved_rows=rr[2][0],sites=len(site_map)))
        # Check explicit annotation/identity failures and role stability with
        # the actually executed helper; all failure outputs stay in this audit.
        unit=args.out/'unit_inputs';unit.mkdir()
        rows=[dict(baseline_examples[name]) for name in ('discovery','reserved')]
        def write(name,values,headers=baseline_fields):
            path=unit/(name+'.tsv')
            with path.open('x',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=headers,delimiter='\t',lineterminator='\n')
                writer.writeheader();writer.writerows(values)
            return path
        checks=[]
        for name,values,headers in (
            ('duplicate_header',rows,baseline_fields+['AD']),
            ('conflicting_alt',[rows[0],dict(rows[0],ALT='G' if rows[0]['ALT']!='G' else 'T')],baseline_fields),
            ('fractional_position',[dict(rows[0],POS='1.5')],baseline_fields),
            ('identical_ref_alt',[dict(rows[0],ALT=rows[0]['REF'])],baseline_fields)):
            path=write(name,values,headers)
            try:
                roles.split_table(path,args.out/(name+'_split'),genome_build='hg38',seed=20261108)
            except ValueError as error:
                checks.append(dict(case=name,error=str(error)))
            else:
                raise AssertionError('Invalid annotation accepted: '+name)
        baseline=roles.split_table(write('unit_baseline',rows),args.out/'unit_baseline_split',genome_build='hg38',seed=20261108)
        changed=[dict(row) for row in rows];changed[1]['GT']='1|0' if changed[1]['GT']=='0|1' else '0|1'
        for name,values in (('reordered',list(reversed(rows))),('reserved_phase',changed)):
            result=roles.split_table(write(name,values),args.out/(name+'_split'),genome_build='hg38',seed=20261108)
            assert result['sites']==baseline['sites']
            assert (args.out/(name+'_split/discovery.tsv')).read_bytes()==(args.out/'unit_baseline_split/discovery.tsv').read_bytes()
        comparisons=[]
        for model in ('binomial','beta'):
            root=args.task/(model+'_v1');base=root/cases[0]
            report=json.loads((root/'report.json').read_text())
            scientific=[p.name for p in base.iterdir() if p.name.endswith(('.tsv','.tsv.gz')) or p.name=='confirmation.json']
            posterior=np.load(base/'posteriors.npy')
            for case in cases:
                directory=root/case
                child=json.loads((directory/'run_manifest.json').read_text())
                assert child['status']=='complete' and child['backend']=='cpu'
                for group in ('input','output'):
                    for record in child[group]['files'].values():
                        assert file_record(record['path'])['sha256']==record['sha256']
                for name in scientific:
                    def contents(path):
                        return gzip.open(path,'rb').read() if path.suffix=='.gz' else path.read_bytes()
                    data=contents(directory/name)
                    assert data==contents(base/name)
                    assert hashlib.sha256(data).hexdigest()==report['fingerprints'][case][name]
                assert np.array_equal(np.load(directory/'posteriors.npy'),posterior,equal_nan=True)
                assert file_record(directory/'posteriors.npy')['sha256']==report['fingerprints'][case]['posteriors.npy']
            cells=pd.read_csv(base/'cells.tsv.gz',sep='\t')
            query=cells[cells.role.eq('tumour_epithelium')&cells.status.eq('assigned')]
            comparisons.append(dict(model=model,raw_runs=4,scientific_tables=len(scientific),
                assigned_query_cells=len(query),MAP_group_sizes=query.groupby('clone').size().to_dict(),
                inference='Candidate/reference-like labels, not validated DNA clones'))
        write_json(dict(conservation=conservation,metadata_rejections=checks,
                        independent_whole_site_hashes_verified=True,role_reorder_phase_append_invariant=True,
                        public_runs_checked=8,models=comparisons,biological_validation=False,
                        error_power_claim=False,production_changed=False),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
