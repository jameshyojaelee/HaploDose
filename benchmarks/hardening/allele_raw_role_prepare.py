"""Prepare locked raw-input role-isolation fixtures, without fitting outcomes."""
import argparse
import importlib.util
import json
from pathlib import Path

import pandas as pd

from haplodose.provenance import start_run,finish_run
from haplodose.simulate import write_dataset


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    helper=Path(__file__).with_name('allele_site_roles.py')
    spec=importlib.util.spec_from_file_location('site_roles',helper);roles=importlib.util.module_from_spec(spec);spec.loader.exec_module(roles)
    run=start_run(args.out,dict(script=Path(__file__),helper=helper,design=args.task/'DESIGN.md'),
                  dict(seed=20261108,query_cells=300,reference_cells_per_type=120),
                  dict(scope='Synthetic raw-input computational fixtures only'),mode='allele_raw_role_prepare')
    try:
        specimens=[dict(name='A',n_cells=300,normal_frac=.2,clones=[
            dict(frac=.5,events=[('chr1',0.,1.,3),('chr2',0.,1.,1)]),
            dict(frac=.5,events=[('chr1',0.,1.,3),('chr2',0.,1.,1),('chr3',0.,1.,2,'loh')])])]
        data=write_dataset(args.out/'data',seed=20261108,n_chrom=4,genes_per_chrom=120,
                           n_types=2,ref_cells_per_type=120,specimens=specimens,alleles=True,
                           snp_rate=1.,snp_capture=.15,ref_patient='A',extra_genes=())
        raw=pd.read_csv(data/'alleles.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
        keys=raw[['CHROM','POS']].drop_duplicates()
        lookup={(r.CHROM,r.POS):roles.site_role(r.CHROM,r.POS,genome_build='hg38',seed=20261108)
                for r in keys.itertuples(index=False)}
        mask=[lookup[(r.CHROM,r.POS)]=='reserved' for r in raw[['CHROM','POS']].itertuples(index=False)]
        variants={'baseline':raw.copy(),'reserved_flip':raw.copy(),'reserved_zero':raw.copy()}
        variants['reserved_flip'].loc[mask,'AD']=(raw.loc[mask,'DP'].astype(int)-raw.loc[mask,'AD'].astype(int)).astype(str)
        variants['reserved_zero'].loc[mask,['AD','DP']]='0'
        positions=pd.read_csv(data/'positions.tsv',sep='\t')
        occupied={roles.canonical_site(r.CHROM,r.POS) for r in keys.itertuples(index=False)}
        extra=None
        for gene in positions.itertuples(index=False):
            for pos in range(int(gene.start),int(gene.end)+1):
                if (gene.chromosome,pos) not in occupied and roles.site_role(gene.chromosome,pos,genome_build='hg38',seed=20261108)=='reserved':
                    extra=raw.iloc[0].copy();extra['CHROM']=gene.chromosome;extra['POS']=str(pos)
                    extra['cM']=str(pos/1e6);extra['gene']=gene.gene;extra['REF']='A';extra['ALT']='C'
                    extra['AD']='1';extra['DP']='2';extra['GT']='0|1';extra['snp_id']='added_metadata_fixture'
                    break
            if extra is not None:
                break
        assert extra is not None
        variants['reserved_added']=pd.concat([raw,pd.DataFrame([extra])],ignore_index=True)
        reports={}
        for name,frame in variants.items():
            path=args.out/(name+'.tsv');frame.to_csv(path,sep='\t',index=False)
            reports[name]=roles.split_table(path,args.out/(name+'_roles'),genome_build='hg38',seed=20261108)
        original=(args.out/'baseline_roles/discovery.tsv').read_bytes()
        assert all((args.out/(name+'_roles/discovery.tsv')).read_bytes()==original for name in variants)
        before=reports['baseline']['sites']
        assert all(reports[name]['sites'][key]==value for name in reports for key,value in before.items())
        # Independently ordered lookup calls and chromosome aliases cannot
        # change existing assignments, without balancing from observed counts.
        for (chrom,pos),role in reversed(list(lookup.items())):
            assert roles.site_role(chrom.removeprefix('chr'),pos,genome_build='hg38',seed=20261108)==role
        (args.out/'report.json').write_text(json.dumps(dict(cases=list(variants),
            unique_original_sites=len(before),reserved_added_site=f"{extra['CHROM']}:{extra['POS']}",
            discovery_files_identical=True,role_order_alias_append_invariant=True,
            biological_specimens=0),indent=2))
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
