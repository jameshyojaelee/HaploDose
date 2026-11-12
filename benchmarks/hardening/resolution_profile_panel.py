"""Fixed physical-resolution RNA profiles on the existing raw-bulk regression."""
import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from haplodose.bins import chromosome_starts, arm_starts
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record
from haplodose.smoothing import predictive_rna_profile, _running_median


def load(name):
    spec = importlib.util.spec_from_file_location('resolution_'+name, Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def grid(start, end, blocks, size=10000000):
    start, end, blocks = np.asarray(start), np.asarray(end), np.asarray(blocks)
    if start.shape != end.shape or start.shape != blocks.shape or start.ndim != 1 or not len(start):
        raise ValueError('Aligned nonempty intervals/blocks required')
    if start.dtype.kind not in 'iu' or end.dtype.kind not in 'iu' or np.any(end <= start):
        raise ValueError('Ordered integer interval endpoints required')
    physical = (start+(end-start)//2)//size
    if np.any((np.diff(start) < 0) & (blocks[1:] == blocks[:-1])):
        raise ValueError('Intervals must be ordered within genomic blocks')
    cut = np.r_[True, (blocks[1:] != blocks[:-1]) | (physical[1:] != physical[:-1])]
    fine_group = np.cumsum(cut)-1
    first = np.flatnonzero(cut)
    # Missing physical intervals terminate smoothers. A declared existing gap
    # also starts a new block even when it lies inside the same physical tile.
    coarse_reset = np.r_[True, (blocks[first[1:]] != blocks[first[:-1]]) |
        (physical[first[1:]] != physical[first[:-1]]+1)]
    coarse_block = np.cumsum(coarse_reset)-1
    return fine_group, coarse_block, first, physical[first]


def checks():
    group, blocks, first, physical = grid(np.array([0,8,10,30,31,0]),
        np.array([2,10,12,31,33,2]), np.array([0,0,0,0,1,2]), size=10)
    np.testing.assert_array_equal(group, [0,0,1,2,3,4])
    np.testing.assert_array_equal(blocks, [0,0,1,2,3])
    np.testing.assert_array_equal(physical, [0,1,3,3,0])
    g = np.repeat(group, 3)
    values = np.repeat(np.array([0.,0.,1.,3.,4.,8.]), 3)
    predicted, record = predictive_rna_profile(values, g, blocks)
    expected = np.array([0.,1.,3.,4.,8.])
    np.testing.assert_array_equal(predicted, expected)
    return dict(independent_boundary_gap_fixture=True, noiseless_profile=True,
        declared_resolution_bp=10000000)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    root = args.task.parents[1]
    native = root/'dev/native_predictive_profile_20261109/raw_bulk_v1/haplodose'
    prior = root/'dev/bulk_amplitude_reference_20261109/refit_000_v1'
    panel = root/'dev/hardening_20261106_2119/profile_v6'
    base = load('profile')
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        frozen=prior/'window_scores.npz', mapping=prior/'gene_map.tsv', split=panel/'split_manifest.tsv',
        scorer=Path(base.__file__), transitive_scorer=Path(base.subarm.__file__),
        statistics=Path(__file__).with_name('statistics.py'))
    for name in ('gene_rna_log_ratio.npy','rna_log_ratio_predictive.npy','genes.tsv.gz','bins.tsv','samples.tsv'):
        inputs[name] = native/name
    run = start_run(args.out, inputs, dict(physical_resolution_bp=10000000),
        dict(scope='Existing raw-bulk regression; fixed RNA grid; no independent evaluation'),
        mode='resolution_profile')
    try:
        numerical = checks()
        bins = pd.read_csv(inputs['bins.tsv'], sep='\t')
        genes = pd.read_csv(inputs['genes.tsv.gz'], sep='\t')
        samples = pd.read_csv(inputs['samples.tsv'], sep='\t', dtype={'sample':str})
        values = np.load(inputs['gene_rna_log_ratio.npy'])
        old = np.load(inputs['rna_log_ratio_predictive.npy'])
        fine_blocks = np.cumsum(chromosome_starts(bins)|arm_starts(bins, genome_build='hg38'))-1
        group, blocks, first, tiles = grid(bins.start.to_numpy(), bins.end.to_numpy(), fine_blocks)
        gene_groups, gene_bins = group[genes.bin.to_numpy()], genes.bin.to_numpy()
        predicted = {name:[] for name in ('physical_10mb_median','physical_10mb_predictive','gene_median_w1','gene_median_w5')}
        records = []
        for j, sample in enumerate(samples['sample']):
            baseline, metadata = predictive_rna_profile(values[j], gene_bins, fine_blocks)
            np.testing.assert_array_equal(baseline, old[j])
            fitted, decision = predictive_rna_profile(values[j], gene_groups, blocks)
            physical = np.array([np.median(values[j, gene_groups == b]) for b in range(len(blocks))])
            medians = np.array([np.median(values[j, gene_bins == b]) for b in range(len(bins))])
            predicted['physical_10mb_median'].append(physical[group])
            predicted['physical_10mb_predictive'].append(fitted[group])
            predicted['gene_median_w1'].append(medians)
            predicted['gene_median_w5'].append(_running_median(medians, fine_blocks, 5))
            records.append(dict(sample=sample, **decision))
        predicted = {k:np.asarray(v) for k,v in predicted.items()}
        predicted['gene_predictive'] = old
        np.savez_compressed(args.out/'RNA_predictions.npz', **predicted, fine_group=group,
            physical_tiles=tiles, first_bin=first, coarse_blocks=blocks)
        write_json(dict(cases=records, predictions=file_record(args.out/'RNA_predictions.npz')),
            args.out/'RNA_DECISIONS_LOCKED.json')
        # Truth opens only after every RNA prediction and parameter is saved.
        with np.load(inputs['frozen']) as z:
            source = {k:z[k] for k in z.files}
        mapping = pd.read_csv(inputs['mapping'], sep='\t')
        originals = [k for k in source if k not in ('samples','common','loss','neutral','gain','dna','split')]
        assert len(originals) == 12
        windows = {name:source[name] for name in originals}
        indexed = genes.set_index('gene').bin
        for name, array in predicted.items():
            aligned = base.subarm.align_rows(array, samples['sample'], source['samples'])
            mapped = base.subarm.gene_values(aligned, indexed, mapping.gene.to_numpy())
            windows[name] = base.subarm.window_median(np.where(mapping.shared.to_numpy(bool)[:,None], mapped, np.nan),
                mapping.window.to_numpy(int), source['dna'].shape[0])
            assert np.isfinite(windows[name][source['common']]).all(), name
        split = pd.read_csv(inputs['split'],sep='\t',dtype=str).set_index('sample').loc[source['samples']]
        assert not split.patient.duplicated().any() and np.array_equal(split.split, source['split'])
        calibration = np.flatnonzero(split.split.eq('calibration'))
        assert len(calibration) == 24
        scores = []
        for name, array in windows.items():
            threshold = base.loss_threshold(array, source['neutral']&source['common'], calibration)
            losses = base.score_loss(array, source['loss'], source['neutral'], source['common'], threshold,
                negative=source['gain']|source['neutral'])
            for j, sample in enumerate(source['samples']):
                x, y = array[source['common'][:,j],j], source['dna'][source['common'][:,j],j]
                cor = float(pearsonr(x,y).statistic) if len(x)>=3 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                if name.startswith('physical_') and len(x)>=3 and np.ptp(y)>0 and not np.isfinite(cor):
                    raise ArithmeticError('Candidate loses an informative patient')
                scores.append(dict(sample=str(sample),split=split.iloc[j].split,method=name,pearson=cor,**losses[j]))
        scores = pd.DataFrame(scores)
        scores.to_csv(args.out/'patient_scores.tsv',sep='\t',index=False)
        scores.groupby(['split','method'])[['pearson','recall_loss','precision_loss','fpr_loss']].agg(['mean','count']).to_csv(args.out/'summary.tsv',sep='\t')
        paired = []
        evaluation = scores[scores.split.eq('evaluation')]
        for metric in ('pearson','recall_loss','precision_loss','fpr_loss'):
            matrix = evaluation.pivot(index='sample',columns='method',values=metric)
            for candidate in ('physical_10mb_median','physical_10mb_predictive'):
                for name in windows:
                    if name == candidate:
                        continue
                    delta = (matrix[candidate]-matrix[name]).dropna()
                    improvement = -delta if metric == 'fpr_loss' else delta
                    paired.append(dict(candidate=candidate,metric=metric,baseline=name,patients=len(delta),
                        mean=float(delta.mean()),wins=int((improvement>0).sum()),losses=int((improvement<0).sum()),
                        ties=int((delta==0).sum()),higher_is_better=metric!='fpr_loss'))
        write_json(dict(checks=numerical,paired=paired,original_comparators=originals,
            fine_bins=len(bins),physical_bins=len(blocks),all_common_support_preserved=True,
            independent_evaluation=False,production_changed=False,default_promoted=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU physical-resolution RNA selection and complete fixed bulk scoring')
    except BaseException as error:
        finish_run(args.out,run,error=error)
        raise


if __name__ == '__main__':
    main()
