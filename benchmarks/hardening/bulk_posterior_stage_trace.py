"""Localize raw bulk signal changes before and after HMM segment averaging."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from haplodose.bulk import run_bulk
from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True); parser.add_argument('--task',type=Path,required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    source = Path(__file__).parent
    spec = importlib.util.spec_from_file_location('fixed_bulk_profile',source/'profile.py')
    base = importlib.util.module_from_spec(spec); spec.loader.exec_module(base)
    panel = args.task.parent/'hardening_20261106_2119/profile_v6'; inputs_dir = panel/'raw_refit_000'
    inputs = dict(script=Path(__file__),design=args.task/'STAGE_TRACE_DESIGN.md',
        profile=source/'profile.py',scorer=source.parent/'tcga_bulk/subarm.py',statistics=source/'statistics.py',
        support=panel/'raw_support_and_truth.npz',genes=panel/'gene_universe.txt',split=panel/'split_manifest.tsv',
        dna_gtf=base.subarm.GTF37,cytoband=base.subarm.CYTO_HG19,centromeres=base.subarm.CEN38,
        absolute_calls=base.subarm.CN/'TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
        absolute_segments=base.subarm.CN/'TCGA_mastercalls.abs_segtabs.fixed.txt')
    for name in ('counts.tsv','samples.tsv','positions.tsv'): inputs[name] = inputs_dir/name
    run = start_run(args.out,inputs,dict(reference_refit=0,min_cpm=0.,scientific_parameter_selection=False),
        dict(scope='Same exposed raw-STAR development panel; descriptive stage attribution only'),mode='bulk_posterior_stage_trace')
    try:
        hd = run_bulk(inputs_dir/'counts.tsv',inputs_dir/'samples.tsv',inputs_dir/'positions.tsv',args.out/'haplodose',
            min_cpm=0.,device='cpu',counts_kind='raw',genome_build='hg38',coordinate_system='one-based-closed',gene_namespace='symbol')
        sample = pd.read_csv(hd/'samples.tsv',sep='\t',dtype={'sample':str},float_precision='round_trip')
        effect = np.load(hd/'rna_log_effect.npy'); state = np.load(hd/'rna_state_index.npy')
        probability = np.load(hd/'rna_state_posterior.npy')
        segments = pd.read_csv(hd/'segments.tsv',sep='\t',dtype={'sample':str},float_precision='round_trip')
        projected = np.full_like(effect,np.nan); coverage = np.zeros(effect.shape,np.int64)
        for row in segments.itertuples():
            i = sample['sample'].tolist().index(row.sample); a,b = row.first_bin,row.last_bin
            assert row.rna_log_effect==float(effect[i,a:b+1].mean())
            assert np.all(state[i,a:b+1]==row.cn_state-1)
            projected[i,a:b+1] = row.rna_log_effect; coverage[i,a:b+1] += 1
        assert (coverage==1).all() and np.isfinite(effect).all() and np.isfinite(probability).all()
        np.save(args.out/'segment_projected_effect.npy',projected)
        # Every RNA fit and stage projection is complete before truth values enter scoring.
        genes = (panel/'gene_universe.txt').read_text().splitlines()
        g37 = base.subarm.hg19_genes(); arms = base.subarm.hg19_arms(); tab,seg = base.subarm.absolute()
        d = base.profiles(inputs_dir,hd,(g37,tab,seg,*arms),genes)
        with np.load(panel/'raw_support_and_truth.npz',allow_pickle=True) as archive:
            frozen = {key:archive[key] for key in archive.files}
        assert np.array_equal(d['tum'],frozen['samples'])
        for key,old in (('D','dna'),('FL','loss'),('NEU','neutral')):
            assert np.array_equal(d[key],frozen[old],equal_nan=True),key
        common = np.isfinite(d['D']) & np.logical_and.reduce([np.isfinite(d['W'][name]) for name in base.METHODS])
        assert np.array_equal(common,frozen['common'])
        mapping = pd.read_csv(hd/'genes.tsv.gz',sep='\t').set_index('gene').bin
        mapped = base.subarm.gene_values(base.subarm.align_rows(effect,sample['sample'],d['tum']),mapping,d['genes'])
        mapped = np.where(d['shared_gene_mapping'][:,None],mapped,np.nan)
        posterior_windows = base.subarm.window_median(mapped,d['window_id'],len(d['warm']))
        assert np.array_equal(common & np.isfinite(posterior_windows),common)
        d['W']['hd_posterior'] = posterior_windows
        # Independent scalar sorting checks every new mapped window, not just favourable loci.
        scalar_checks = 0
        for w in range(len(d['warm'])):
            for j in range(len(d['tum'])):
                values = sorted(float(x) for x in mapped[d['window_id']==w,j] if np.isfinite(x))
                if not values: assert np.isnan(posterior_windows[w,j]); continue
                n = len(values); expected = values[n//2] if n%2 else (values[n//2-1]+values[n//2])/2
                assert posterior_windows[w,j]==np.float32(expected)
                scalar_checks += 1
        split = pd.read_csv(panel/'split_manifest.tsv',sep='\t',dtype=str).set_index('sample').loc[d['tum']]
        ci = np.flatnonzero(split.split.eq('calibration')); ei = np.flatnonzero(split.split.eq('evaluation'))
        assert len(ci)==len(ei)==24 and not split.patient.duplicated().any()
        rows = []
        for method in (*base.METHODS,'hd_posterior'):
            values = d['W'][method]
            threshold = base.loss_threshold(values,d['NEU']&common,ci)
            loss = base.score_loss(values,d['FL'],d['NEU'],common,threshold,negative=d['FG']|d['NEU'])
            for j,name in enumerate(d['tum']):
                mask = common[:,j]; x,y = values[mask,j],d['D'][mask,j]
                valid = len(x)>=3 and np.ptp(x)>0 and np.ptp(y)>0
                focal = common[:,j]&d['FL'][:,j]
                rows.append(dict(sample=str(name),split=str(split.iloc[j].split),method=method,
                    callable_windows=int(mask.sum()),rho=float(spearmanr(x,y).statistic) if valid else np.nan,
                    pearson=float(pearsonr(x,y).statistic) if valid else np.nan,
                    loss_windows=int(focal.sum()),loss_signal_median=float(np.median(values[focal,j])) if focal.any() else np.nan,
                    loss_threshold=threshold,**loss[j]))
        scores = pd.DataFrame(rows); scores.to_csv(args.out/'patient_stage_scores.tsv',sep='\t',index=False)
        np.savez_compressed(args.out/'window_stages.npz',samples=d['tum'].astype(str),common=common,
            posterior=posterior_windows,segment=d['W']['hd_dose'],ratio=d['W']['hd_lr'],loss=d['FL'],neutral=d['NEU'])
        summaries = scores.groupby(['split','method']).agg(
            n_samples=('sample','size'),correlation_patients=('pearson','count'),mean_pearson=('pearson','mean'),
            loss_patients=('recall_loss','count'),mean_loss_recall=('recall_loss','mean'),mean_loss_fpr=('fpr_loss','mean'))
        summaries.to_csv(args.out/'descriptive_summary.tsv',sep='\t')
        report = dict(samples=len(sample),tumour_patients=len(d['tum']),reference_samples=int(sample.role.eq('normal').sum()),
            scalar_window_checks=scalar_checks,segments_exact=len(segments),old_callable_support_exact=True,
            original_comparator_family_retained=list(base.METHODS),added_stage='hd_posterior',
            truth='Existing uncertain gene-anchored ABSOLUTE DNA windows; no new biological replication',
            threshold_source='Same24 calibration patients only; no jitter',inferential_tests_added=0,
            stage_summary='descriptive_summary.tsv',method_promotion=False,scientific_criteria_changed=False)
        (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error)
        raise


if __name__=='__main__': main()
