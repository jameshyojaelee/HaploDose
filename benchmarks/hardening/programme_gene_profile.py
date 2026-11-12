"""Research programme-adjusted gene RNA profiles with exact native controls."""
import argparse
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import statistics

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from haplodose import bulk
from haplodose.bins import chromosome_starts, arm_starts
from haplodose.io import write_json
from haplodose.profile import ProfileParams, log_ratio
from haplodose.provenance import start_run, finish_run, file_record
from haplodose.smoothing import predictive_rna_profile, _running_median


def load(name):
    spec = importlib.util.spec_from_file_location('programme_gene_'+name, Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gene_expectation(expected, reference, gene_bins):
    expected, reference, gene_bins = np.asarray(expected, float), np.asarray(reference, float), np.asarray(gene_bins, int)
    totals = np.bincount(gene_bins, weights=reference, minlength=len(expected))
    if (reference.shape != gene_bins.shape or not np.isfinite(reference).all() or np.any(reference < 0)
            or np.any(totals <= 0) or not np.isfinite(expected).all() or np.any(expected <= 0)):
        raise ValueError('Supported reference composition and positive bin expectations required')
    values = expected[gene_bins]*(reference/totals[gene_bins])
    rebuilt = np.bincount(gene_bins, weights=values, minlength=len(expected))
    operations = 4*len(reference)+32
    unit = np.finfo(float).eps
    budget = operations*unit/(1-operations*unit)*max(1.,float(expected.max()))
    if np.max(abs(rebuilt-expected)) > budget:
        raise ArithmeticError('Gene expectations do not conserve bin expectations')
    return values


def region_medians(values, gene_bins, path, reset):
    path, reset = np.asarray(path), np.asarray(reset, bool)
    if len(path) != len(reset) or not reset[0]:
        raise ValueError('Aligned state path with initial reset required')
    cut = reset.copy()
    cut[1:] |= path[1:] != path[:-1]
    group = np.cumsum(cut)-1
    output = np.empty(len(path))
    for k in np.unique(group):
        bins = np.flatnonzero(group == k)
        output[bins] = np.median(values[np.isin(gene_bins,bins)])
    return output


def checks():
    reference = np.array([3.,7.,2.,6.])
    bins = np.array([0,0,1,1])
    expectation = gene_expectation([20.,16.],reference,bins)
    exact = [Fraction(20)*Fraction(3,10),Fraction(20)*Fraction(7,10),
             Fraction(16)*Fraction(2,8),Fraction(16)*Fraction(6,8)]
    np.testing.assert_array_equal(expectation,list(map(float,exact)))
    np.testing.assert_array_equal(log_ratio(expectation,expectation),np.zeros(4))
    values = np.array([0.,8.,2.,4.,10.,12.])
    mapping = np.array([0,0,1,1,2,3])
    fitted = region_medians(values,mapping,[0,0,0,0],[True,False,True,False])
    independent = [statistics.median([0.,8.,2.,4.])]*2+[statistics.median([10.,12.])]*2
    np.testing.assert_array_equal(fitted,independent)
    permutation = np.array([3,0,5,2,1,4])
    np.testing.assert_array_equal(fitted,region_medians(values[permutation],mapping[permutation],
        [0,0,0,0],[True,False,True,False]))
    return dict(independent_fraction_expectations=True, count_unit_identity=True,
        scalar_region_medians=True, gene_permutation_and_gap_resets=True)


def panel(args):
    root = args.task.parents[1]
    number = args.reference_refit
    native = (root/'dev/native_predictive_profile_20261109/raw_bulk_v1/haplodose' if number == 0 else
        root/f'dev/hardening_20261106_2119/profile_v6/raw_refit_{number:03d}/haplodose')
    frozen = root/f'dev/bulk_amplitude_reference_20261109/refit_{number:03d}_v1'
    split_path = root/'dev/hardening_20261106_2119/profile_v6/split_manifest.tsv'
    prior_manifest = json.loads((native/'run_manifest.json').read_text())
    base = load('profile')
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        native_manifest=native/'run_manifest.json', frozen=frozen/'window_scores.npz',
        mapping=frozen/'gene_map.tsv', split=split_path, scorer=Path(base.__file__),
        transitive_scorer=Path(base.subarm.__file__), statistics=Path(__file__).with_name('statistics.py'))
    for key, record in prior_manifest['input']['files'].items():
        inputs[key] = Path(record['path'])
    names = ['log_ratio.npy','genes.tsv.gz','bins.tsv','samples.tsv','run.json']
    if number == 0:
        names += ['gene_rna_log_ratio.npy','rna_log_ratio_predictive.npy','rna_state_index.npy']
    for name in names:
        inputs[name] = native/name
    run = start_run(args.out,inputs,dict(candidate_family=2,reference_refit=number),
        dict(scope='Existing raw-bulk regression; conditional programme adjustment; no independent validation'),
        mode='programme_gene_profile')
    try:
        numerical = checks()
        for record in prior_manifest['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        genes = pd.read_csv(inputs['genes.tsv.gz'],sep='\t')
        bins = pd.read_csv(inputs['bins.tsv'],sep='\t')
        samples = pd.read_csv(inputs['samples'],sep='\t',dtype=str)
        native_samples = pd.read_csv(inputs['samples.tsv'],sep='\t',dtype={'sample':str},float_precision='round_trip')
        np.testing.assert_array_equal(samples['sample'],native_samples['sample'])
        table = pd.read_csv(inputs['counts'],sep='\t',index_col=0)
        assert table.index.is_unique and genes.gene.is_unique
        np.testing.assert_array_equal(table.index.to_numpy()[genes.column.to_numpy(int)],genes.gene)
        G = table.loc[:,samples['sample']].to_numpy(float)[genes.column.to_numpy(int)]
        gbin = genes.bin.to_numpy(int)
        Y = np.stack([G[gbin == b].sum(axis=0) for b in range(len(bins))],axis=1)
        count_column = 'input_binned_counts' if 'input_binned_counts' in native_samples else 'binned_counts'
        np.testing.assert_array_equal(Y.sum(axis=1),native_samples[count_column])
        normal = samples.role.eq('normal').to_numpy()
        reset = chromosome_starts(bins)
        arm_reset = arm_starts(bins,genome_build='hg38')
        blocks = np.cumsum(reset|arm_reset)-1
        parameters = json.loads(inputs['run.json'].read_text())
        params = ProfileParams(**parameters['params'])
        n_programmes = parameters.get('n_programmes_requested',parameters['n_programmes'])
        full_ref = bulk.reference_model(Y[normal],n_programmes,arm_reset)
        old_gene = np.load(inputs['gene_rna_log_ratio.npy']) if number == 0 else None
        old_profile = np.load(inputs['rna_log_ratio_predictive.npy']) if number == 0 else None
        old_lr = np.load(inputs['log_ratio.npy'])
        old_path = np.load(inputs['rna_state_index.npy']) if number == 0 else None
        members = [np.flatnonzero(gbin == b) for b in range(len(bins))]
        totals = G.sum(axis=0)
        ref_gene, ref_depth = G[:,normal].sum(axis=1), float(totals[normal].sum())
        predictions = {k:[] for k in ('programme_gene_predictive','programme_gene_native_regions',
            'raw_gene_native_regions','gene_median_w1','gene_median_w5','gene_predictive')}
        decisions, contexts = [], []
        for i,sample in enumerate(samples['sample']):
            ref = bulk.reference_model(Y[normal & (np.arange(len(Y)) != i)],
                n_programmes,arm_reset) if normal[i] else full_ref
            fit = bulk.fit_sample(Y[i],*ref,reset,params)
            lr = log_ratio(Y[i],fit['E'])
            lr -= np.median(lr)
            np.testing.assert_array_equal(lr.astype(np.float32),old_lr[i])
            if number == 0:
                np.testing.assert_array_equal(fit['path'],old_path[i])
            source_gene = ref_gene-G[:,i] if normal[i] else ref_gene
            source_depth = ref_depth-totals[i] if normal[i] else ref_depth
            raw,centre = bulk._direct_gene_ratio(G[:,i],totals[i]*(source_gene/source_depth),members)
            if number == 0:
                np.testing.assert_array_equal(raw,old_gene[i])
            unchanged,meta = predictive_rna_profile(raw,gbin,blocks)
            if number == 0:
                np.testing.assert_array_equal(unchanged,old_profile[i])
            predictions['gene_predictive'].append(unchanged)
            expected = gene_expectation(fit['E'],source_gene,gbin)
            corrected,new_centre = bulk._direct_gene_ratio(G[:,i],expected,members)
            fitted,decision = predictive_rna_profile(corrected,gbin,blocks)
            predictions['programme_gene_predictive'].append(fitted)
            predictions['programme_gene_native_regions'].append(region_medians(corrected,gbin,fit['path'],reset|arm_reset))
            predictions['raw_gene_native_regions'].append(region_medians(raw,gbin,fit['path'],reset|arm_reset))
            raw_bins = np.array([np.median(raw[m]) for m in members])
            predictions['gene_median_w1'].append(raw_bins)
            predictions['gene_median_w5'].append(_running_median(raw_bins,blocks,5))
            contexts.append(dict(sample=sample,centre=new_centre,programme_shift_rms=float(np.sqrt(np.mean(fit['shift']**2)))))
            np.savez_compressed(args.out/f'context_{i:03d}.npz',expected_gene=expected,
                expected_bin=fit['E'],corrected_gene=corrected,shift=fit['shift'],path=fit['path'])
            decisions.append(dict(sample=sample,**decision))
        predictions = {k:np.asarray(v) for k,v in predictions.items()}
        np.savez_compressed(args.out/'RNA_predictions.npz',**predictions)
        write_json(dict(cases=decisions,contexts=contexts,predictions=file_record(args.out/'RNA_predictions.npz'),
            exact_native_controls=len(samples)),args.out/'RNA_DECISIONS_LOCKED.json')
        # All nuisance, selection and boundaries are saved before truth opens.
        with np.load(inputs['frozen']) as z:
            truth = {k:z[k] for k in z.files}
        mapping = pd.read_csv(inputs['mapping'],sep='\t')
        originals = [k for k in truth if k not in ('samples','common','loss','neutral','gain','dna','split')]
        assert len(originals) == 12
        windows = {k:truth[k] for k in originals}
        for name,array in predictions.items():
            aligned = base.subarm.align_rows(array,samples['sample'],truth['samples'])
            mapped = base.subarm.gene_values(aligned,genes.set_index('gene').bin,mapping.gene.to_numpy())
            windows[name] = base.subarm.window_median(np.where(mapping.shared.to_numpy(bool)[:,None],mapped,np.nan),
                mapping.window.to_numpy(int),truth['dna'].shape[0])
            assert np.isfinite(windows[name][truth['common']]).all(),name
        np.savez_compressed(args.out/'window_scores.npz',**truth,**{k:v for k,v in windows.items() if k not in truth})
        split = pd.read_csv(split_path,sep='\t',dtype=str).set_index('sample').loc[truth['samples']]
        assert not split.patient.duplicated().any() and np.array_equal(split.split,truth['split'])
        calibration = np.flatnonzero(split.split.eq('calibration'))
        assert len(calibration) == 24
        scores = []
        for name,array in windows.items():
            threshold = base.loss_threshold(array,truth['neutral']&truth['common'],calibration)
            loss = base.score_loss(array,truth['loss'],truth['neutral'],truth['common'],threshold,
                negative=truth['gain']|truth['neutral'])
            for j,sample in enumerate(truth['samples']):
                x,y = array[truth['common'][:,j],j],truth['dna'][truth['common'][:,j],j]
                cor = float(pearsonr(x,y).statistic) if len(x)>=3 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                if name.startswith(('programme_gene_','raw_gene_native_')) and len(x)>=3 and np.ptp(y)>0 and not np.isfinite(cor):
                    raise ArithmeticError('Candidate loses an informative patient')
                scores.append(dict(sample=str(sample),split=split.iloc[j].split,method=name,pearson=cor,**loss[j]))
        scores = pd.DataFrame(scores)
        scores.to_csv(args.out/'patient_scores.tsv',sep='\t',index=False)
        scores.groupby(['split','method'])[['pearson','recall_loss','precision_loss','fpr_loss']].agg(['mean','count']).to_csv(args.out/'summary.tsv',sep='\t')
        paired = []
        selected = scores[scores.split.eq('evaluation')]
        for metric in ('pearson','recall_loss','precision_loss','fpr_loss'):
            matrix = selected.pivot(index='sample',columns='method',values=metric)
            for candidate in ('programme_gene_predictive','programme_gene_native_regions'):
                for baseline in windows:
                    if baseline == candidate:
                        continue
                    delta = (matrix[candidate]-matrix[baseline]).dropna()
                    improvement = -delta if metric == 'fpr_loss' else delta
                    paired.append(dict(candidate=candidate,baseline=baseline,metric=metric,patients=len(delta),
                        mean=float(delta.mean()),wins=int((improvement>0).sum()),losses=int((improvement<0).sum()),
                        ties=int((delta==0).sum()),higher_is_better=metric!='fpr_loss'))
        write_json(dict(checks=numerical,exact_native_controls=len(samples),methods=list(windows),paired=paired,
            independent_evaluation=False,production_changed=False,default_promoted=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU native bulk fits and programme-adjusted gene RNA profiles')
    except BaseException as error:
        finish_run(args.out,run,error=error)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--checks-only',action='store_true')
    parser.add_argument('--reference-refit',type=int,choices=range(9),default=0)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    if args.checks_only:
        run = start_run(args.out,dict(script=Path(__file__),design=args.task/'DESIGN.md'),{},
            dict(scope='Independent scalar fixtures'),mode='programme_gene_numerics')
        try:
            write_json(checks(),args.out/'report.json')
            finish_run(args.out,run,backend='CPU scalar/Fraction fixtures')
        except BaseException as error:
            finish_run(args.out,run,error=error)
            raise
    else:
        panel(args)


if __name__ == '__main__':
    main()
