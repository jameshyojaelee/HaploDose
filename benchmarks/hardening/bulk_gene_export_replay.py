"""Require the package gene effect to reproduce all frozen direct-gene RNA inputs."""
import argparse
import hashlib
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.bulk import _direct_gene_ratio
from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    root=Path('/gpfs/commons/home/jameslee/Haplodose')
    panel=root/'dev/hardening_20261106_2119/profile_v6'
    helper=root/'dev/gene_aggregation_20261108/candidate_v2/benchmarks/hardening/gene_aggregation_panel.py'
    spec=importlib.util.spec_from_file_location('frozen_direct_gene',helper)
    reference=importlib.util.module_from_spec(spec);spec.loader.exec_module(reference)
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',reference=helper)
    for i in range(9):
        for name in ('counts.tsv','samples.tsv','haplodose/genes.tsv.gz'):
            inputs[f'{i}:{name}']=panel/f'raw_refit_{i:03d}'/name
    run=start_run(args.out,inputs,{},dict(scope='Exact development gene-effect replay, not new biological evaluation'),
                  mode='bulk_gene_export_replay')
    try:
        records=[]
        for i in range(9):
            directory=panel/f'raw_refit_{i:03d}'
            counts=pd.read_csv(directory/'counts.tsv',sep='\t',index_col=0)
            meta=pd.read_csv(directory/'samples.tsv',sep='\t',dtype=str).set_index('sample').loc[counts.columns]
            mapping=pd.read_csv(directory/'haplodose/genes.tsv.gz',sep='\t')
            assert counts.index.is_unique and mapping.gene.is_unique and counts.columns.is_unique
            G=counts.loc[mapping.gene].to_numpy(float);normal=meta.role.eq('normal').to_numpy()
            bin_ids=mapping.bin.to_numpy(int)
            expected_reference,_=reference.gene_effects(G,bin_ids,normal)
            members=[np.flatnonzero(bin_ids==j) for j in range(bin_ids.max()+1)]
            library=G.sum(axis=0);rate=G[:,normal].sum(axis=1)/library[normal].sum()
            targets=np.flatnonzero(~normal);digest=hashlib.sha256()
            for sample in targets:
                value,_=_direct_gene_ratio(G[:,sample],library[sample]*rate,members)
                np.testing.assert_array_equal(value,expected_reference[:,sample])
                digest.update(value.tobytes())
            records.append(dict(refit=i,tumour_profiles=len(targets),genes=len(G),
                                exact_gene_effect_identity=True,package_effect_sha256=digest.hexdigest()))
            print(f'Exact gene-effect replay {i+1}/9',flush=True)
        write_json(dict(refits=records,profiles=sum(r['tumour_profiles'] for r in records),
                        prospective_evaluation=False,superiority_claim=False),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
