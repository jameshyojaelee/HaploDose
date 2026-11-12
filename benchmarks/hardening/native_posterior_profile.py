"""Require installed bulk RNA kernels to reproduce retained predictions exactly."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.bins import chromosome_starts, arm_starts
from haplodose.io import write_json
from haplodose.posterior import posterior_rna_profile, average_rna_profiles
from haplodose.profile import ProfileParams
from haplodose.provenance import start_run, finish_run, file_record


METADATA=('samples','common','loss','neutral','gain','dna','split')


def central(args,root):
    prior=root/'dev/programme_gene_profile_20261109/refit_000_v1'
    frozen=root/'dev/posterior_region_average_20261109/bulk_v1'
    source=json.loads((prior/'run_manifest.json').read_text())['input']['files']
    inputs=dict(script=Path(__file__),design=args.task/'PORT_DESIGN.md',
        prior=prior/'run_manifest.json',frozen=frozen/'RNA_predictions.npz',
        decisions=frozen/'Gaussian_decisions.json')
    for name in('samples','genes.tsv.gz','bins.tsv','run.json'):
        inputs[name]=Path(source[name]['path'])
    for i in range(68):inputs[f'context{i}']=prior/f'context_{i:03d}.npz'
    rec=start_run(args.out,inputs,dict(scope='All48 central profiles before DNA scoring'),
        dict(calibration='Descriptive RNA only; no new biological validation'),mode='native_posterior_central_parity')
    try:
        genes=pd.read_csv(inputs['genes.tsv.gz'],sep='\t')
        bins=pd.read_csv(inputs['bins.tsv'],sep='\t')
        samples=pd.read_csv(inputs['samples'],sep='\t',dtype=str)
        assert len(samples)==68 and genes.gene.is_unique and samples['sample'].is_unique
        normal=samples.role.eq('normal').to_numpy();assert normal.sum()==20
        effects,expected=[],[]
        for i in range(68):
            with np.load(inputs[f'context{i}'])as z:
                effects.append(z['corrected_gene']);expected.append(z['expected_gene'])
        effects,expected=map(np.asarray,(effects,expected))
        assert effects.shape==expected.shape==(68,len(genes))
        reset=chromosome_starts(bins)|arm_starts(bins,genome_build='hg38')
        blocks=np.cumsum(reset)-1
        params=ProfileParams(**json.loads(inputs['run.json'].read_text())['params'])
        decisions={r['sample']:r for r in json.loads(inputs['decisions'].read_text())}
        with np.load(inputs['frozen'])as z:old=z['averaged_gaussian_reference_median']
        output=np.full(old.shape,np.nan);records=[]
        for i in np.flatnonzero(~normal):
            value,decision=posterior_rna_profile(effects[i],expected[i],genes.bin.to_numpy(int),
                effects[normal],expected[normal],blocks,reference_id='frozen-central-20-normal-panel',
                counts_kind='raw',params=params)
            np.testing.assert_array_equal(value,old[i])
            sample=str(samples.iloc[i]['sample']);before=decisions[sample]
            for key in('beta','tau','centre','grid_scores','near_tied_grid_indices','score_rounding_budget'):
                assert decision['boundary'][key]==before[key],(sample,key)
            output[i]=value;records.append(dict(sample=sample,row=int(i),exact=True,decision=decision))
            write_json(dict(sample=sample,exact=True),args.out/f'case_{i:03d}.json')
            print(json.dumps(dict(sample=sample,profile_exact=True)),flush=True)
        np.testing.assert_array_equal(output,old)
        np.savez_compressed(args.out/'RNA_predictions.npz',profiles=output,samples=samples['sample'].to_numpy(str))
        write_json(records,args.out/'decisions.json')
        write_json(dict(artifact=file_record(args.out/'RNA_predictions.npz'),DNA_accessed=False,
            complete_profiles_exact=48,reference_rows_unchanged=True),args.out/'RNA_LOCKED.json')
        write_json(dict(complete_profiles_exact=48,boundary_decisions_exact=48,
            DNA_accessed=False,scientific_model_changed=False,defaults_changed=False,
            biological_validation=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU installed public posterior RNA profile, complete48 central exact replay')
    except BaseException as error:
        finish_run(args.out,rec,error=error);raise


def ensemble(args,root):
    inputs=dict(script=Path(__file__),design=args.task/'PORT_DESIGN.md',
        expected=root/'dev/reference_bagged_profile_20261110/aggregate_v1/RNA_windows.npz',
        expected_lock=root/'dev/reference_bagged_profile_20261110/aggregate_v1/RNA_LOCKED.json')
    for i in range(1,9):inputs[f'refit{i}']=root/f'dev/posterior_region_average_20261109/refit_{i:03d}_v1/window_scores.npz'
    rec=start_run(args.out,inputs,dict(refits=list(range(1,9)),methods=39),
        dict(scope='All fixed RNA window means; no DNA archive members accessed'),mode='native_posterior_ensemble_parity')
    try:
        assert json.loads(inputs['expected_lock'].read_text())['DNA_arrays_accessed']is False
        fits,identifiers={},None
        for i in range(1,9):
            with np.load(inputs[f'refit{i}'])as z:
                names=[n for n in z.files if n not in METADATA]
                if identifiers is None:
                    assert len(names)==39
                    identifiers=z['samples'];fits={name:[]for name in names}
                np.testing.assert_array_equal(identifiers,z['samples']);assert names==list(fits)
                for name in names:fits[name].append(z[name])
        outputs={}
        with np.load(inputs['expected'])as old:
            np.testing.assert_array_equal(identifiers,old['samples'])
            for name,arrays in fits.items():
                windows,patients=arrays[0].shape
                assert all(a.shape==(windows,patients)for a in arrays)
                axis=np.arange(windows);ids=[f'refit_{i:03d}'for i in range(1,9)]
                output=np.full_like(arrays[0],np.nan,dtype=float)
                for j in range(patients):
                    output[:,j],_=average_rna_profiles([a[:,j]for a in arrays],
                        [axis for _ in arrays],reference_ids=ids)
                key='reference_bag_'+name
                np.testing.assert_array_equal(output,old[key]);outputs[key]=output
        np.savez_compressed(args.out/'RNA_windows.npz',samples=identifiers,**outputs)
        write_json(dict(artifact=file_record(args.out/'RNA_windows.npz'),methods_exact=39,
            DNA_arrays_accessed=False,operator='Identical equal mean AFTER mapping each fit to common windows'),
            args.out/'RNA_LOCKED.json')
        write_json(dict(methods_exact=39,reference_fits=8,complete_specimen_columns=len(identifiers),
            DNA_arrays_accessed=False,scientific_model_changed=False,biological_validation=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU installed identity-aligned reference ensemble, all39 exact RNA arrays')
    except BaseException as error:
        finish_run(args.out,rec,error=error);raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--mode',choices=('central','ensemble'),required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    (central if args.mode=='central'else ensemble)(args,Path(os.environ['HAPLODOSE_DATA_ROOT']))


if __name__=='__main__':main()
