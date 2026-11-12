"""Run stable-SNP-role raw-input isolation through the unmodified public pipeline."""
import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path

import jax

from haplodose.alleles import AlleleParams
from haplodose.clones import CloneParams
from haplodose.io import write_json
from haplodose.pipeline import run as run_pipeline
from haplodose.provenance import start_run,finish_run,file_record


def outputs(directory):
    result={}
    for path in sorted(directory.iterdir()):
        if path.name.endswith(('.tsv','.tsv.gz')) or path.name=='confirmation.json':
            data=gzip.open(path,'rb').read() if path.suffix=='.gz' else path.read_bytes()
            result[path.name]=hashlib.sha256(data).hexdigest()
        elif path.name=='posteriors.npy':
            result[path.name]=file_record(path)['sha256']
    qc=copy.deepcopy(json.loads((directory/'qc.json').read_text()))
    for name in ('started','wall_seconds','max_rss_gib'):
        qc.pop(name)
    qc['alleles'].pop('path')
    result['numerical_qc']=hashlib.sha256(json.dumps(qc,sort_keys=True).encode()).hexdigest()
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--cell-model',choices=('binomial','beta'),required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    assert jax.default_backend()=='cpu'
    source=args.task/'prepare_v1'
    manifest=json.loads((source/'run_manifest.json').read_text())
    assert manifest['status']=='complete'
    for group in ('input','output'):
        for record in manifest[group]['files'].values():
            assert file_record(record['path'])['sha256']==record['sha256']
    report=json.loads((source/'report.json').read_text())
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',prepare=source/'run_manifest.json')
    for name in ('counts.h5','cells.tsv','positions.tsv'):
        inputs[name]=source/'data'/name
    for case in report['cases']:
        inputs[case]=source/(case+'_roles/discovery.tsv')
    record=start_run(args.out,inputs,dict(cell_model=args.cell_model,cases=report['cases']),
                     dict(scope='Raw public-run computational isolation; no biological/coverage/power claim'),
                     mode='allele_raw_role_run')
    try:
        fingerprints={};summary=[]
        for case in report['cases']:
            print(json.dumps(dict(case=case,model=args.cell_model,status='starting')),flush=True)
            destination=args.out/case
            run_pipeline(inputs['counts.h5'],inputs['cells.tsv'],inputs['positions.tsv'],destination,
                         alleles=inputs[case],params=CloneParams(seed=20261108,split_mode='thin'),
                         allele_params=AlleleParams(cell_kappa=None if args.cell_model=='binomial' else 2.),
                         counts_kind='raw',genome_build='hg38',coordinate_system='one-based-closed',
                         gene_namespace='symbol',allele_genome_build='hg38',device='cpu')
            child=json.loads((destination/'run_manifest.json').read_text())
            assert child['status']=='complete' and child['backend']=='cpu'
            fingerprint=outputs(destination)
            fingerprints[case]=fingerprint
            assert fingerprint==fingerprints[report['cases'][0]], 'Reserved-only raw input changed scientific outputs'
            qc=json.loads((destination/'qc.json').read_text())
            summary.append(dict(case=case,n_bins=qc['n_bins'],n_genes=qc['n_genes_used'],
                                snps=qc['alleles']['n_snps'],scientific_outputs_identical=True))
        write_json(dict(models=args.cell_model,cases=summary,fingerprints=fingerprints,
                        scope='Research role splitter plus public raw-input discovery; reserved observations not fitted',
                        biological_specimens=0,calibrated=False,production_changed=False),args.out/'report.json')
        for item in record['input']['files'].values():
            assert file_record(item['path'])['sha256']==item['sha256']
        finish_run(args.out,record,backend='cpu')
    except Exception as error:
        finish_run(args.out,record,error=error);raise


if __name__=='__main__':
    main()
