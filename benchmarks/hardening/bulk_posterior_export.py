"""Exact raw-count/CLI delivery gate for optional bulk posterior RNA profiles."""
import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose import bulk, cli
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


ARRAYS=('log_ratio.npy','gene_rna_log_ratio.npy','rna_state_posterior.npy',
        'rna_log_effect.npy','rna_state_index.npy')
TABLES=('bins.tsv','genes.tsv.gz','segments.tsv','arms.tsv','samples.tsv','gene_ratio_metadata.tsv')


def exact_native(before,after):
    for name in ARRAYS:
        np.testing.assert_array_equal(np.load(before/name),np.load(after/name))
    for name in TABLES:
        pd.testing.assert_frame_equal(pd.read_csv(before/name,sep='\t',float_precision='round_trip'),
            pd.read_csv(after/name,sep='\t',float_precision='round_trip'),check_exact=True)
    return dict(arrays=len(ARRAYS),tables=len(TABLES),exact=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--mode',choices=('control','candidate'),required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    root=Path(os.environ['HAPLODOSE_DATA_ROOT'])
    frozen=root/'dev/native_predictive_profile_20261109/raw_bulk_v1/haplodose'
    prior=root/'dev/programme_gene_profile_20261109/refit_000_v1'
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',
        old_manifest=frozen/'run_manifest.json',old_run=frozen/'run.json')
    source=json.loads(inputs['old_manifest'].read_text())
    config=json.loads(inputs['old_run'].read_text())
    assert config['backend']=='cpu' and config['bin_stat']=='sum'
    assert config['input_contract']['counts_kind']=='raw'
    assert json.loads(json.dumps(asdict(bulk.BULK_PROFILE)))==config['params']
    for name in('counts','samples','positions'):
        inputs[name]=Path(source['input']['files'][name]['path'])
        assert file_record(inputs[name])==source['input']['files'][name]
    for name in ARRAYS+TABLES:inputs['old_'+name]=frozen/name
    if args.mode=='candidate':
        inputs['control']=args.task/'control_v1/report.json'
        inputs['control_manifest']=args.task/'control_v1/run_manifest.json'
        assert json.loads(inputs['control_manifest'].read_text())['status']=='complete'
        assert json.loads(inputs['control'].read_text())['frozen_native_exact']['exact']
        inputs['frozen_RNA']=root/'dev/posterior_region_average_20261109/bulk_v1/RNA_predictions.npz'
        for i in range(68):inputs[f'context{i}']=prior/f'context_{i:03d}.npz'
    out=args.out/'bulk'
    argv=['bulk','--counts',str(inputs['counts']),'--samples',str(inputs['samples']),
        '--positions',str(inputs['positions']),'--out',str(out),'--device','cpu',
        '--bin-genes',str(config['bin_genes']),'--min-cpm',str(config['min_cpm']),
        '--programmes',str(config.get('n_programmes_requested',config['n_programmes'])),
        '--bin-stat','sum','--rna-profile-smoothing','none'if args.mode=='control'else'posterior-regions']
    for key in('counts_kind','genome_build','coordinate_system','gene_namespace'):
        argv+=['--'+key.replace('_','-'),str(config['input_contract'][key])]
    if config['input_contract'].get('centromeres')is not None:
        argv+=['--centromeres',config['input_contract']['centromeres']]
    if config.get('max_gap_bp')is not None:argv+=['--max-gap-bp',str(config['max_gap_bp'])]
    if 'exclude_genes'in source['input']['files']:
        argv+=['--exclude-genes',source['input']['files']['exclude_genes']['path']]
    rec=start_run(args.out,inputs,dict(mode=args.mode,cli_argv=argv),
        dict(scope='Complete raw-count/CLI RNA delivery; no DNA access'),mode='bulk_posterior_export')
    try:
        # Exercise the installed argument parser and the actual public bulk
        # entry point, not a research replacement of fit_sample or exports.
        cli.main(argv)
        assert json.loads((out/'run_manifest.json').read_text())['status']=='complete'
        assert json.loads((out/'run.json').read_text())['backend']=='cpu'
        native=exact_native(frozen,out)
        report=dict(frozen_native_exact=native,actual_inference_backend='cpu',
            mode=args.mode,DNA_accessed=False,scientific_model_changed=False,defaults_changed=False)
        if args.mode=='candidate':
            report['control_native_exact']=exact_native(args.task/'control_v1/bulk',out)
            expected=np.load(out/'gene_fitted_expected.npy')
            effect=np.load(out/'gene_rna_log_ratio_programme.npy')
            samples=pd.read_csv(out/'samples.tsv',sep='\t',dtype=str)
            assert len(samples)==68
            for i in range(68):
                with np.load(inputs[f'context{i}'])as z:
                    np.testing.assert_array_equal(expected[i],z['expected_gene'])
                    np.testing.assert_array_equal(effect[i],z['corrected_gene'])
            normal=samples.role.eq('normal').to_numpy();assert normal.sum()==20
            value=np.load(out/'rna_log_ratio_posterior_regions.npy')
            with np.load(inputs['frozen_RNA'])as z:
                np.testing.assert_array_equal(value,z['averaged_gaussian_reference_median'])
            assert np.isnan(value[normal]).all()and np.isfinite(value[~normal]).all()
            write_json(dict(artifact=file_record(out/'rna_log_ratio_posterior_regions.npy'),
                complete_tumour_profiles_exact=48,reference_only_rows=20,
                complete_gene_contexts_exact=68,DNA_accessed=False),args.out/'RNA_LOCKED.json')
            report.update(complete_tumour_profiles_exact=48,complete_gene_contexts_exact=68,
                reference_only_rows=20,smoother_backend='cpu',calibration_status='unvalidated')
        write_json(report,args.out/'report.json')
        finish_run(args.out,rec,backend='CPU installed full raw-count bulk CLI')
    except BaseException as error:
        finish_run(args.out,rec,error=error);raise


if __name__=='__main__':main()
