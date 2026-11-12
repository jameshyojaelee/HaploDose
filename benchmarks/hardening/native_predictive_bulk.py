"""Native raw-bulk candidate versus the complete frozen comparator family."""
import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from haplodose.bins import chromosome_starts, arm_starts
from haplodose.bulk import run_bulk
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record
from haplodose.smoothing import _running_median


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    root = Path('/gpfs/commons/home/jameslee/Haplodose')
    panel = root/'dev/hardening_20261106_2119/profile_v6'
    inp = panel/'raw_refit_000'
    prior = root/'dev/bulk_posterior_export_20261109/stage_v1/haplodose'
    frozen = root/'dev/bulk_amplitude_reference_20261109/refit_000_v1/window_scores.npz'
    reference = root/'dev/bulk_amplitude_reference_20261109/refit_000_v1/gene_map.tsv'
    scorer = Path(__file__).with_name('profile.py')
    spec = importlib.util.spec_from_file_location('predictive_bulk_scoring', scorer)
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    inputs = dict(script=Path(__file__), design=args.task/'RAW_BULK_DESIGN.md',
        scorer=scorer, transitive_scorer=Path(base.subarm.__file__),
        statistics=Path(__file__).with_name('statistics.py'), frozen=frozen, mapping=reference,
        dna_gtf=base.subarm.GTF37, cytoband=base.subarm.CYTO_HG19, centromeres=base.subarm.CEN38,
        absolute_calls=base.subarm.CN/'TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
        absolute_segments=base.subarm.CN/'TCGA_mastercalls.abs_segtabs.fixed.txt',
        split=panel/'split_manifest.tsv', universe=panel/'gene_universe.txt')
    for name in ('counts.tsv', 'samples.tsv', 'positions.tsv'):
        inputs[name] = inp/name
    for path in prior.iterdir():
        if path.suffix == '.npy' or path.name.endswith(('.tsv', '.tsv.gz')):
            inputs['frozen_native:'+path.name] = path
    run = start_run(args.out, inputs, dict(default_change=False, reference_refit=0),
        dict(scope='Former evaluation/development regression only; no new biological replication'),
        mode='native_predictive_bulk')
    try:
        out = run_bulk(inp/'counts.tsv', inp/'samples.tsv', inp/'positions.tsv', args.out/'haplodose',
            min_cpm=0., device='cpu', counts_kind='raw', genome_build='hg38',
            coordinate_system='one-based-closed', gene_namespace='symbol', rna_profile_smoothing='predictive')
        for name, path in inputs.items():
            if not name.startswith('frozen_native:'):
                continue
            if path.suffix == '.npy':
                np.testing.assert_array_equal(np.load(path), np.load(out/path.name))
            else:
                pd.testing.assert_frame_equal(pd.read_csv(path,sep='\t'),
                                              pd.read_csv(out/path.name,sep='\t'),check_exact=True)
        bins = pd.read_csv(out/'bins.tsv', sep='\t')
        genes = pd.read_csv(out/'genes.tsv.gz', sep='\t')
        samples = pd.read_csv(out/'samples.tsv', sep='\t', dtype={'sample': str})
        blocks = np.cumsum(chromosome_starts(bins) | arm_starts(bins, genome_build='hg38'))-1
        direct = np.load(out/'gene_rna_log_ratio.npy')
        medians = np.array([[np.median(row[genes.bin == b]) for b in range(len(bins))] for row in direct])
        methods = dict(gene_median_w1=medians,
            gene_median_w5=np.array([_running_median(row,blocks,5) for row in medians]),
            gene_predictive=np.load(out/'rna_log_ratio_predictive.npy'))
        # All candidate RNA values are now fixed, before outcome values are loaded.
        write_json({name:file_record(out/name) for name in
            ('rna_log_ratio_predictive.npy','rna_profile_smoothing.json','gene_rna_log_ratio.npy')},
            args.out/'RNA_PREDICTIONS_LOCKED.json')
        truth = base.subarm.hg19_genes()
        calls, segments = base.subarm.absolute()
        d = base.profiles(inp,out,(truth,calls,segments,*base.subarm.hg19_arms()),
                          (panel/'gene_universe.txt').read_text().splitlines())
        with np.load(frozen) as z:
            old = {key:z[key] for key in z.files}
        np.testing.assert_array_equal(d['tum'],old['samples'])
        common = np.isfinite(d['D']) & np.logical_and.reduce([np.isfinite(d['W'][m]) for m in base.METHODS])
        for value, key in ((common,'common'),(d['D'],'dna'),(d['FL'],'loss'),(d['NEU'],'neutral'),(d['FG'],'gain')):
            np.testing.assert_array_equal(value,old[key])
        names = [name for name in old if name not in ('samples','common','loss','neutral','gain','dna','split')]
        assert len(names) == 12
        windows = {name:old[name] for name in names}
        mapping = genes.set_index('gene').bin
        for name,values in methods.items():
            mapped = base.subarm.gene_values(base.subarm.align_rows(values,samples['sample'],d['tum']),mapping,d['genes'])
            windows[name] = base.subarm.window_median(np.where(d['shared_gene_mapping'][:,None],mapped,np.nan),d['window_id'],len(d['warm']))
            assert np.isfinite(windows[name][common]).all(), name
        split = pd.read_csv(panel/'split_manifest.tsv',sep='\t',dtype=str).set_index('sample').loc[d['tum']]
        assert not split.patient.duplicated().any()
        calibration = np.flatnonzero(split.split.eq('calibration'))
        assert len(calibration) == 24 and int(split.split.eq('evaluation').sum()) == 24
        rows = []
        for name,values in windows.items():
            threshold = base.loss_threshold(values,d['NEU']&common,calibration)
            losses = base.score_loss(values,d['FL'],d['NEU'],common,threshold,negative=d['FG']|d['NEU'])
            for j,sample in enumerate(d['tum']):
                x,y = values[common[:,j],j], d['D'][common[:,j],j]
                correlation = float(pearsonr(x,y).statistic) if len(x)>=3 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                rows.append(dict(sample=str(sample), split=split.iloc[j].split, method=name,
                                 pearson=correlation, **losses[j]))
        scores = pd.DataFrame(rows)
        scores.to_csv(args.out/'patient_scores.tsv',sep='\t',index=False)
        scores.groupby(['split','method'])[['pearson','recall_loss','precision_loss','fpr_loss',
            'classified_truth_coverage']].agg(['mean','count']).to_csv(args.out/'descriptive_summary.tsv',sep='\t')
        paired = []
        evaluation = scores[scores.split.eq('evaluation')]
        for metric in ('pearson','recall_loss','precision_loss','fpr_loss'):
            matrix = evaluation.pivot(index='sample',columns='method',values=metric)
            for baseline in [*names,'gene_median_w1','gene_median_w5']:
                delta = (matrix.gene_predictive-matrix[baseline]).dropna()
                paired.append(dict(metric=metric,baseline=baseline,patients=len(delta),
                    mean_difference=float(delta.mean()),wins=int((delta>0).sum()),
                    losses=int((delta<0).sum()),ties=int((delta==0).sum()),
                    higher_is_better=metric!='fpr_loss'))
        write_json(dict(native_arrays_and_tables_unchanged=True, paired_descriptive=paired,
            original_comparators=names, independent_validation=False, default_promoted=False,
            scientific_criteria_changed=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU native raw-count bulk and resolution-matched frozen DNA regression')
    except BaseException as error:
        finish_run(args.out,run,error=error)
        raise


if __name__ == '__main__':
    main()
