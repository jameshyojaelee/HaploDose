"""Trace every callable comparison loss window to counts and local RNA evidence.

The independent-bin posterior is a diagnostic counterfactual, not a new caller.
No parameters, thresholds, truth labels or support rules are selected here.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import softmax

from haplodose.bins import arm_starts, chromosome_starts
from haplodose.bulk import BULK_PROFILE, fit_sample, reference_model
from haplodose.device import select
from haplodose.kernels.likelihood import pln_logpmf
from haplodose.profile import LOG_DOSE, STATE_PRIOR, log_ratio
from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    root = args.task.parent.parent
    prior = root/'dev/bulk_posterior_export_20261109/stage_v1'
    hd = prior/'haplodose'
    panel = root/'dev/hardening_20261106_2119/profile_v6'
    inp = panel/'raw_refit_000'
    source = Path(__file__).parent
    spec = importlib.util.spec_from_file_location('fixed_bulk_trace_profile', source/'profile.py')
    base = importlib.util.module_from_spec(spec); spec.loader.exec_module(base)
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        profile=source/'profile.py', scorer=source.parent/'tcga_bulk/subarm.py',
        split=panel/'split_manifest.tsv', genes=panel/'gene_universe.txt',
        dna_gtf=base.subarm.GTF37, cytoband=base.subarm.CYTO_HG19, centromeres=base.subarm.CEN38,
        absolute_calls=base.subarm.CN/'TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
        absolute_segments=base.subarm.CN/'TCGA_mastercalls.abs_segtabs.fixed.txt')
    for directory, prefix, names in ((inp, 'input', ('counts.tsv','samples.tsv','positions.tsv')),
        (prior, 'stage', ('window_stages.npz','patient_stage_scores.tsv','run_manifest.json')),
        (hd, 'hd', ('bins.tsv','genes.tsv.gz','samples.tsv','log_ratio.npy',
                    'rna_state_posterior.npy','rna_log_effect.npy','rna_state_index.npy'))):
        inputs.update({prefix+':'+name:directory/name for name in names})
    run = start_run(args.out, inputs, dict(parameter_selection=False, truth_selection=False),
        dict(scope='Previously exposed panel; all callable comparison loss windows; descriptive mechanism trace'),
        mode='bulk_likelihood_trace')
    try:
        backend = select('cpu'); assert backend=='cpu'
        bins = pd.read_csv(hd/'bins.tsv', sep='\t', float_precision='round_trip')
        gene = pd.read_csv(hd/'genes.tsv.gz', sep='\t')
        meta = pd.read_csv(hd/'samples.tsv', sep='\t', float_precision='round_trip')
        counts = pd.read_csv(inp/'counts.tsv', sep='\t', index_col=0).loc[gene.gene, meta['sample']]
        assert counts.index.is_unique and np.equal(counts, np.floor(counts)).all().all()
        aggregation = sparse.csr_matrix((np.ones(len(gene)), (np.arange(len(gene)), gene.bin)),
            shape=(len(gene),len(bins)))
        Y = np.asarray((aggregation.T @ counts.to_numpy(np.float64)).T)
        assert np.array_equal(Y.sum(axis=1), meta.input_binned_counts)
        normal = meta.role.eq('normal').to_numpy()
        reset = chromosome_starts(bins)
        mu, v, W = reference_model(Y[normal], 10, arm_starts(bins, genome_build='hg38'))
        assert np.array_equal(np.sqrt(v),bins.reference_sd.to_numpy())
        with np.load(prior/'window_stages.npz') as z: stage = {k:z[k] for k in z.files}
        split = pd.read_csv(panel/'split_manifest.tsv', sep='\t').set_index('sample').loc[stage['samples']]
        selected = stage['loss'] & stage['common'] & split.split.eq('evaluation').to_numpy()[None]
        assert selected.sum()==49
        jselected = np.flatnonzero(selected.any(axis=0))
        scores = pd.read_csv(prior/'patient_stage_scores.tsv', sep='\t', float_precision='round_trip')
        thresholds = {m:float(scores.loc[scores.method.eq(m),'loss_threshold'].unique().item())
            for m in ('hd_lr','hd_posterior')}
        g37 = base.subarm.hg19_genes(); tab, seg = base.subarm.absolute()
        d = base.profiles(inp,hd,(g37,tab,seg,*base.subarm.hg19_arms()),
            (panel/'gene_universe.txt').read_text().splitlines())
        assert np.array_equal(d['tum'],stage['samples'])
        for key,old in (('FL','loss'),('NEU','neutral')): assert np.array_equal(d[key],stage[old])
        assert np.array_equal(d['W']['hd_lr'],stage['ratio'],equal_nan=True)
        mapping = gene.set_index('gene').bin.reindex(d['genes'])
        post = np.load(hd/'rna_state_posterior.npy'); dose = np.load(hd/'rna_log_effect.npy')
        path = np.load(hd/'rna_state_index.npy'); ratio = np.load(hd/'log_ratio.npy')
        rows = []; gene_rows = []; numerical = []; fits = {}
        for j in jselected:
            sample = str(stage['samples'][j]); i = meta['sample'].tolist().index(sample)
            f = fit_sample(Y[i],mu,v,W,reset,BULK_PROFILE)
            for key,expected in (('post',post[i]),('dose',dose[i]),('path',path[i])):
                assert np.array_equal(f[key],expected), (sample,key)
            for key in ('centre','beta','tau'): assert f[key]==meta.iloc[i][key]
            r = log_ratio(Y[i],f['E']); median = np.median(r)
            assert np.array_equal((r-median).astype(np.float32),ratio[i])
            means = np.log(f['E'].astype(np.float32))[:,None] + np.float32(f['centre']) + f['beta']*LOG_DOSE
            emit = np.asarray(pln_logpmf(Y[i].astype(np.float32)[:,None],means,f['sd'].astype(np.float32)[:,None]))
            local = softmax(emit.astype(float)+np.log(STATE_PRIOR.astype(float)),axis=1)
            fits[sample] = dict(Y=Y[i],E=f['E'],shift=f['shift'],sd=f['sd'],emit=emit,
                means=means,local=local,centre=f['centre'],beta=f['beta'],tau=f['tau'],
                raw_model_residual=r-f['centre'],median_centre=median)
            np.savez_compressed(args.out/(sample+'.npz'),**fits[sample])
            for w in np.flatnonzero(selected[:,j]):
                gi = np.flatnonzero((d['window_id']==w)&d['shared_gene_mapping'])
                bi = mapping.iloc[gi].to_numpy(int); unique = np.unique(bi)
                raw_call = stage['ratio'][w,j]<thresholds['hd_lr']
                posterior_call = stage['posterior'][w,j]<thresholds['hd_posterior']
                category = 'both' if raw_call and posterior_call else 'raw_only' if raw_call else 'posterior_only' if posterior_call else 'neither'
                # All medians weight genes exactly as the fixed window scorer.
                values = dict(sample=sample,window=int(w),category=category,arm=str(d['warm'][w]),
                    hg38_start=int(d['gene_positions'].iloc[gi].start.min()),
                    hg38_last_gene_start=int(d['gene_positions'].iloc[gi].start.max()),
                    genes=len(gi),unique_bins=len(unique),dna_relative_log2=float(d['D'][w,j]),
                    raw_score=float(stage['ratio'][w,j]),posterior_score=float(stage['posterior'][w,j]),
                    local_effect=float(np.float32(np.median((local@(f['beta']*LOG_DOSE))[bi]))),
                    local_loss_probability=float(np.median(local[bi,0])),
                    smoothed_loss_probability=float(np.median(post[i,bi,0])),
                    local_loss_neutral_log_bayes_factor=float(np.median(emit[bi,0]-emit[bi,1])),
                    context_loss_neutral_log_odds_change=float(np.median(
                        np.log(np.maximum(post[i,bi,0],np.finfo(np.float32).tiny))-
                        np.log(np.maximum(post[i,bi,1],np.finfo(np.float32).tiny))-
                        (emit[bi,0]-emit[bi,1]+np.log(float(STATE_PRIOR[0]/STATE_PRIOR[1]))))),
                    median_raw_model_residual=float(np.median((r-f['centre'])[bi])),
                    median_reference_sd=float(np.median(np.sqrt(v[bi]))),
                    median_total_sd=float(np.median(f['sd'][bi])),
                    median_programme_shift=float(np.median(f['shift'][bi])),
                    centre=f['centre'],median_centre=median,beta=f['beta'],tau=f['tau'],
                    loss_state_effect=float(f['beta']*LOG_DOSE[0]),
                    raw_count_unique_bins=float(Y[i,unique].sum()),
                    expected_count_unique_bins=float(f['E'][unique].sum()))
                rows.append(values)
                for g,b in zip(gi,bi): gene_rows.append(dict(sample=sample,window=int(w),
                    gene=str(d['genes'][g]),bin=int(b),chromosome=str(bins.iloc[b].chromosome),
                    bin_start=int(bins.iloc[b].start),bin_end=int(bins.iloc[b].end)))
            print(json.dumps(dict(sample=sample,fit_exact=True,loss_windows=int(selected[:,j].sum()))),flush=True)
        frame = pd.DataFrame(rows); frame.to_csv(args.out/'loss_windows.tsv',sep='\t',index=False)
        pd.DataFrame(gene_rows).to_csv(args.out/'window_gene_bins.tsv',sep='\t',index=False)
        frame.groupby('category').agg(windows=('window','size'),patients=('sample','nunique'),
            median_local_loss=('local_loss_probability','median'),median_smoothed_loss=('smoothed_loss_probability','median'),
            median_loss_neutral_log_bf=('local_loss_neutral_log_bayes_factor','median'),
            median_local_effect=('local_effect','median'),median_posterior_effect=('posterior_score','median'),
            median_beta=('beta','median'),median_tau=('tau','median')).to_csv(args.out/'category_summary.tsv',sep='\t')
        report = dict(loss_windows=len(frame),patients=len(jselected),exact_refits=len(fits),
            categories=frame.category.value_counts().to_dict(),thresholds=thresholds,
            reference_donors=int(normal.sum()),biological_replication_added=0,
            local_scope='Independent-bin RNA posterior using identical emissions and state prior; not a calibrated caller',
            DNA_scope='Uncertain ABSOLUTE GRCh37 gene-anchored windows; RNA hg38 positions, not interval liftover',
            parameter_changes=0,method_promotion=False,calibration_claim=False)
        (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        for record in run['input']['files'].values(): assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(args.out,run,backend=backend)
    except Exception as error:
        finish_run(args.out,run,error=error); raise


if __name__=='__main__': main()
