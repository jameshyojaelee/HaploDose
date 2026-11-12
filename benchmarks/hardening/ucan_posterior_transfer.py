"""Challenge the fixed bulk posterior estimator on two locked U-CAN cases."""
import argparse
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from haplodose import bulk, cli
from haplodose.bins import arm_labels
from haplodose.io import write_json
from haplodose.profile import ProfileParams, log_ratio
from haplodose.provenance import start_run, finish_run, file_record


def load(name):
    path = Path(__file__).with_name(name+'.py')
    spec = importlib.util.spec_from_file_location('ucan_transfer_'+name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    original = root/'dev/hardening_20261106_2119/ucan_v6'
    ucan, operator = load('ucan'), load('profile_operator_control')
    provenance = json.loads((original/'provenance.json').read_text())
    assert provenance['locked_samples'] == list(ucan.LOCKED_SAMPLES)
    assert provenance['source_validation']['n_normal_donors'] == 120
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
                  old_provenance=original/'provenance.json',
                  ucan=Path(ucan.__file__), operator=Path(operator.__file__),
                  arm_scoring=Path(__file__).with_name('arm_metrics.py'))
    for name in ('counts.tsv', 'samples.tsv', 'positions.tsv'):
        inputs[name] = original/'input'/name
    for i, record in enumerate(provenance['source_files']):
        inputs[f'raw_source_{i}'] = Path(record['path'])
        assert file_record(record['path']) == record
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, inputs, dict(locked_samples=list(ucan.LOCKED_SAMPLES),
                    smoothing_width_bins=5, programmes=10, bin_genes=10, min_cpm=1.),
                    dict(scope='Locked development transfer, not independent evaluation',
                         biological_query_specimens=2, normal_donors=120),
                    mode='ucan_posterior_transfer')
    try:
        operator.verify()
        X, metadata, positions, verified = ucan.verify_inputs(
            root/'results/ucan_bulk_input_20261101_1352', ucan.RAW_SOURCE, ucan.ANNOTATION)
        X, positions, identity = ucan.canonicalize_observations(X, positions)
        pd.testing.assert_frame_equal(X, pd.read_csv(inputs['counts.tsv'], sep='\t', index_col=0), check_exact=True)
        pd.testing.assert_frame_equal(metadata, pd.read_csv(inputs['samples.tsv'], sep='\t', dtype=str), check_exact=True)
        pd.testing.assert_frame_equal(positions, pd.read_csv(inputs['positions.tsv'], sep='\t', dtype={'chromosome':str}), check_exact=True)
        outputs = {}
        for mode in ('none', 'posterior-regions'):
            out = args.out/mode
            argv = ['bulk', '--counts', str(inputs['counts.tsv']), '--samples', str(inputs['samples.tsv']),
                    '--positions', str(inputs['positions.tsv']), '--out', str(out), '--device', 'cpu',
                    '--bin-genes', '10', '--min-cpm', '1', '--programmes', '10', '--bin-stat', 'sum',
                    '--exclude-genes', str(ucan.MASK), '--counts-kind', 'raw', '--genome-build', 'hg38',
                    '--coordinate-system', 'one-based-closed', '--gene-namespace', 'symbol',
                    '--rna-profile-smoothing', mode]
            write_json(dict(argv=argv, DNA_accessed=False), args.out/(mode+'_argv.json'))
            cli.main(argv)
            assert json.loads((out/'run_manifest.json').read_text())['status'] == 'complete'
            assert json.loads((out/'run.json').read_text())['backend'] == 'cpu'
            outputs[mode] = out
        control, candidate = outputs['none'], outputs['posterior-regions']
        array_names = ('log_ratio.npy', 'gene_rna_log_ratio.npy', 'rna_state_posterior.npy',
                       'rna_log_effect.npy', 'rna_state_index.npy')
        table_names = ('bins.tsv', 'genes.tsv.gz', 'segments.tsv', 'arms.tsv', 'samples.tsv', 'gene_ratio_metadata.tsv')
        for name in array_names:
            np.testing.assert_array_equal(np.load(control/name), np.load(candidate/name))
        for name in table_names:
            pd.testing.assert_frame_equal(pd.read_csv(control/name, sep='\t', float_precision='round_trip'),
                pd.read_csv(candidate/name, sep='\t', float_precision='round_trip'), check_exact=True)
        bins = pd.read_csv(candidate/'bins.tsv', sep='\t')
        used = pd.read_csv(candidate/'genes.tsv.gz', sep='\t')
        samples = pd.read_csv(candidate/'samples.tsv', sep='\t', dtype=str)
        np.testing.assert_array_equal(samples['sample'], metadata['sample'])
        normal = samples.role.eq('normal').to_numpy()
        posterior = np.load(candidate/'rna_log_ratio_posterior_regions.npy')
        assert normal.sum() == 120 and np.isnan(posterior[normal]).all()
        assert np.isfinite(posterior[~normal]).all()
        gene = X.loc[used.gene, samples['sample']].to_numpy(float)
        gene_bin = used.bin.to_numpy(int)
        full_totals = X.sum(axis=0).to_numpy(float)
        logged = np.log(gene/full_totals[None]*1e6+1)
        gene_mu = gene[:,normal].sum(axis=1)/gene[:,normal].sum()
        arms = np.asarray(arm_labels(bins))
        blocks = bins.genomic_block.to_numpy()
        selected, locked = {}, []
        segments = pd.read_csv(candidate/'segments.tsv', sep='\t')
        for sample in ucan.LOCKED_SAMPLES:
            i = int(np.flatnonzero(samples['sample'].eq(sample))[0])
            naive = logged[:,i]-logged[:,normal].mean(axis=1)
            naive -= np.median(naive)
            ratio = log_ratio(gene[:,i],gene[:,i].sum()*gene_mu)
            ratio -= np.median(ratio)
            exported = np.full(len(bins),np.nan)
            for row in segments[segments['sample'].eq(sample)].itertuples():
                exported[row.first_bin:row.last_bin+1] = row.rna_log_effect
            methods = dict(naive=ucan.reduce_genes(naive,gene_bin,len(bins)),
                gene_ratio=ucan.reduce_genes(ratio,gene_bin,len(bins)),
                hd_ratio=np.load(candidate/'log_ratio.npy')[i],
                hd_posterior=np.load(candidate/'rna_log_effect.npy')[i],
                hd_segment_export=exported, posterior_regions=posterior[i])
            assert all(np.isfinite(v).all() for v in methods.values())
            differences = np.concatenate([np.diff(methods['naive'][blocks == b])
                                          for b in pd.unique(blocks) if np.sum(blocks == b)>1])
            sigma = np.median(abs(differences-np.median(differences)))/(np.sqrt(2)*norm.ppf(.75))
            assert np.isfinite(sigma) and sigma > 0
            penalty = 3*sigma**2*np.log(len(bins))
            for name, values in list(methods.items()):
                methods[name+'_smoothed'] = ucan.smooth(values,blocks,5)
                fitted = np.full(len(bins),np.nan)
                for block in pd.unique(blocks):
                    take = blocks == block
                    fitted[take] = operator.partition(values[take],penalty,ProfileParams().min_bins)[0]
                methods[name+'_own_partition'] = fitted
            path = args.out/(sample+'_RNA.npz')
            np.savez_compressed(path,**methods)
            selected[sample] = methods
            locked.append(dict(sample=sample, RNA=file_record(path), sigma=float(sigma),
                               penalty=float(penalty), methods=list(methods)))
        # Both specimens' entire RNA family locks before either DNA file opens.
        write_json(dict(profiles=locked,DNA_accessed=False),args.out/'RNA_LOCKED.json')
        scores, dna_inputs = [], []
        for sample, methods in selected.items():
            patient = metadata.loc[metadata['sample'].eq(sample),'patient'].item()
            dna, header, source = ucan.read_dna(patient)
            dna_inputs.append(file_record(source))
            truth = ucan.reduce_genes(ucan.dna_at_genes(used,dna,header['ploidy']),gene_bin,len(bins),require_complete=True)
            common = np.isfinite(truth)&np.logical_and.reduce([np.isfinite(v) for v in methods.values()])
            frame = bins.copy();frame['truth']=truth;frame['common_callable']=common
            for name, values in methods.items():
                frame[name] = values
                scores.append(dict(sample=sample,method=name,resolution='common_bin',**ucan.metrics(values[common],truth[common])))
                for aggregation in ('median','mean'):
                    scores.append(dict(sample=sample,method=name,resolution='common_arm_'+aggregation,
                        **ucan.arm_metrics(values[common],truth[common],arms[common],aggregation)))
            frame.to_csv(args.out/(sample+'_profiles.tsv.gz'),sep='\t',index=False)
        pd.DataFrame(scores).to_csv(args.out/'scores.tsv',sep='\t',index=False)
        write_json(dict(source_validation=verified,native_arrays_exact=len(array_names),
            native_tables_exact=len(table_names),scores=scores,DNA_sources=dna_inputs,
            biological_query_specimens=2,normal_donors=120,independent_evaluation=False,
            scientific_criteria_changed=False,default_promotion=False),args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out,run,backend='CPU installed raw-count bulk CLI, fixed posterior and simple controls')
    except BaseException as error:
        finish_run(args.out,run,error=error)
        raise


if __name__ == '__main__':
    main()
