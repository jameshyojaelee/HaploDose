"""Test loss-state mass as a ranking score on frozen raw bulk RNA fits."""
import argparse
from fractions import Fraction
import importlib.util
import json
from pathlib import Path

import jax
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from haplodose.io import write_json
from haplodose.provenance import file_record, finish_run, start_run


METHOD = 'hd_posterior_loss_mass'
BASELINES = ('hd_posterior', 'averaged_gaussian_reference_median', 'naive_own_partition')
META = ('samples', 'common', 'loss', 'neutral', 'gain', 'dna', 'split')


def scalar_median(values):
    ordered = sorted(Fraction(float(value)) for value in values if np.isfinite(value))
    n = len(ordered)
    return np.nan if n == 0 else np.float32(float((ordered[(n-1)//2]+ordered[n//2])/2))


def make_rna(args, root, base):
    subarm = base.subarm
    panel = root/'dev/hardening_20261106_2119/profile_v6'
    stage = root/'dev/bulk_posterior_export_20261109/stage_v1'
    hd = stage/'haplodose'
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md', scorer=Path(base.__file__),
        subarm=Path(subarm.__file__), statistics=Path(base.adjust_family.__code__.co_filename),
        posterior=hd/'rna_state_posterior.npy', effect=hd/'rna_log_effect.npy',
        fit_manifest=hd/'run.json', fit_samples=hd/'samples.tsv', fit_genes=hd/'genes.tsv.gz',
        fit_bins=hd/'bins.tsv', stage=stage/'window_stages.npz',
        count_gene_order=panel/'raw_refit_000/counts.tsv', positions=panel/'raw_refit_000/positions.tsv',
        genes=panel/'gene_universe.txt', coordinate_gtf=subarm.GTF37, centromeres=subarm.CEN38)
    rec = start_run(args.out, inputs, dict(primary=METHOD, aggregation='fixed gene-weighted window median'),
                    dict(scope='RNA score and frozen coordinate data; no DNA outcomes decoded'),
                    mode='bulk_loss_state_rna')
    try:
        meta = json.loads(inputs['fit_manifest'].read_text())
        assert meta['bin_rna_state_outputs']['state_labels'] == ['loss','neutral','gain','amplification']
        pos = pd.read_csv(inputs['positions'], sep='\t').set_index('gene')
        pos.index = subarm.canonical_symbols(pos.index, 'RNA position symbols')
        pos['chromosome'] = 'chr'+pos.chromosome.astype(str).str.removeprefix('chr')
        g37 = subarm.hg19_genes()
        g37.index = subarm.canonical_symbols(g37.index, 'DNA position symbols')
        ids = pd.read_csv(inputs['count_gene_order'], sep='\t', usecols=[0]).iloc[:,0]
        ids = pd.Index(subarm.canonical_symbols(ids, 'count genes'))
        fixed = subarm.canonical_symbols(inputs['genes'].read_text().splitlines(), 'fixed scoring symbols')
        keep = ids.isin(fixed) & ids.isin(pos.index) & ids.isin(g37.index)
        keep &= g37.reindex(ids).chromosome.eq(pos.chromosome.reindex(ids)).to_numpy()
        positions = pos.loc[ids[keep]].copy()
        positions['chrom_n'] = positions.chromosome.str[3:].astype(int)
        positions = positions.sort_values(['chrom_n', 'start'])
        cen = pd.read_csv(inputs['centromeres'], sep='\t', comment='#').set_index('chromosome').p_arm_end
        positions['arm'] = positions.chromosome + np.where(positions.start >= positions.chromosome.map(cen), 'q', 'p')
        wid = subarm.windows(positions)
        n_windows = int(wid.max()+1)
        samples = pd.read_csv(inputs['fit_samples'], sep='\t')['sample']
        bins = pd.read_csv(inputs['fit_bins'], sep='\t')['bin']
        np.testing.assert_array_equal(bins.to_numpy(), np.arange(len(bins)))
        mapping = pd.read_csv(inputs['fit_genes'], sep='\t').set_index('gene').bin
        mapping.index = subarm.canonical_symbols(mapping.index, 'HaploDose mapped symbols')
        post = np.load(inputs['posterior'])
        effect = np.load(inputs['effect'])
        assert post.shape == (*effect.shape, 4) and effect.shape == (len(samples), len(bins))
        assert np.isfinite(post).all() and np.all((post >= 0) & (post <= 1))
        with np.load(inputs['stage']) as archive:
            tum, old_effect = archive['samples'], archive['posterior']
        def mapped(x):
            return subarm.gene_values(subarm.align_rows(x, samples, tum), mapping, positions.index)
        np.testing.assert_array_equal(subarm.window_median(mapped(effect), wid, n_windows), old_effect)
        gene_score = mapped(-post[..., 0])
        windows = subarm.window_median(gene_score, wid, n_windows)
        for w in range(n_windows):
            reference = [scalar_median(column) for column in gene_score[wid == w].T]
            np.testing.assert_array_equal(windows[w], reference)
        np.savez_compressed(args.out/'RNA_windows.npz', samples=tum, score=windows,
                            genes=positions.index.to_numpy(dtype=str), window_id=wid)
        write_json(dict(primary=METHOD, identity=file_record(args.out/'RNA_windows.npz'),
            old_posterior_effect_projection_exact=True, scalar_windows_checked=int(windows.size),
            DNA_outcomes_accessed=False), args.out/'RNA_LOCKED.json')
        for item in rec['input']['files'].values():
            assert file_record(item['path'])['sha256']==item['sha256']
        finish_run(args.out, rec, backend='CPU fixed RNA posterior loss-mass score')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


def score(args, root, base, audit):
    previous = root/'dev/robust_reference_profile_20261110/score_v1'
    panel = root/'dev/hardening_20261106_2119/profile_v6'
    rna = args.out.parent/'rna'
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        profile=Path(base.__file__), subarm=Path(base.subarm.__file__),
        statistics=Path(base.adjust_family.__code__.co_filename), reference=Path(audit.__file__),
        windows=previous/'window_scores.npz', old_scores=previous/'patient_scores.tsv',
        RNA=rna/'RNA_windows.npz', RNA_lock=rna/'RNA_LOCKED.json',
        split=panel/'split_manifest.tsv')
    rec = start_run(args.out, inputs, dict(primary=METHOD, fixed_comparators=BASELINES),
        dict(scope='Exposed 48-patient regression; original calibration and evaluation roles'),
        mode='bulk_loss_state_score')
    try:
        lock = json.loads(inputs['RNA_lock'].read_text())
        assert lock['DNA_outcomes_accessed'] is False and lock['identity']==file_record(inputs['RNA'])
        with np.load(inputs['windows']) as z:
            truth = {k:z[k] for k in z.files}
        old = [key for key in truth if key not in META]
        assert len(old)==123 and all(key in old for key in BASELINES)
        with np.load(inputs['RNA']) as z:
            np.testing.assert_array_equal(z['samples'], truth['samples'])
            loss_score = z['score']
        common = truth['common']
        assert loss_score.shape==common.shape and np.isfinite(loss_score[common]).all()
        split = pd.read_csv(inputs['split'], sep='\t', dtype=str).set_index('sample').loc[truth['samples']]
        assert not split.patient.duplicated().any()
        np.testing.assert_array_equal(split.split, truth['split'])
        calibration = np.flatnonzero(split.split.eq('calibration'))
        assert len(calibration)==len(split)-len(calibration)==24
        methods = {key:truth[key] for key in old}
        methods[METHOD] = loss_score
        rows = []
        for name, values in methods.items():
            threshold = base.loss_threshold(values, truth['neutral']&common, calibration)
            losses = base.score_loss(values, truth['loss'], truth['neutral'], common, threshold,
                                     negative=truth['gain']|truth['neutral'])
            if name==METHOD:
                assert threshold==audit.threshold_reference(values, truth['neutral']&common, calibration)
            profile = truth['hd_posterior'] if name==METHOD else values
            for j, sample in enumerate(truth['samples']):
                x, y = profile[common[:,j],j], truth['dna'][common[:,j],j]
                corr = float(pearsonr(x,y).statistic) if len(x)>=3 and np.isfinite(x).all() and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                if name==METHOD:
                    for endpoint, mask in [('recall_loss',truth['loss']), ('fpr_loss',truth['neutral'])]:
                        observed = [float(v) for v in values[common[:,j]&mask[:,j],j]]
                        expected = sum(v < threshold for v in observed)/len(observed) if observed else np.nan
                        actual = losses[j][endpoint]
                        assert actual==expected or (np.isnan(actual) and np.isnan(expected))
                rows.append(dict(sample=str(sample), split=split.iloc[j].split, method=name,
                                 pearson=corr, **losses[j]))
        frame = pd.DataFrame(rows)
        frame.to_csv(args.out/'patient_scores.tsv', sep='\t', index=False)
        replay = pd.read_csv(args.out/'patient_scores.tsv', sep='\t', float_precision='round_trip')
        previous_scores = pd.read_csv(inputs['old_scores'], sep='\t', float_precision='round_trip')
        pd.testing.assert_frame_equal(previous_scores,replay[replay.method.isin(old)].reset_index(drop=True),check_exact=True)
        metrics = ['pearson','recall_loss','fpr_loss','precision_loss']
        evaluated = frame[frame.split.eq('evaluation')]
        summary = frame.groupby(['split','method'])[metrics].agg(['mean','count'])
        summary.to_csv(args.out/'summary.tsv',sep='\t')
        contrasts = {}
        for name in old:
            a = evaluated[evaluated.method.eq(METHOD)].set_index('sample')
            b = evaluated[evaluated.method.eq(name)].set_index('sample')
            np.testing.assert_array_equal(a.index,b.index)
            contrasts[name] = {}
            for endpoint in metrics:
                delta = (a[endpoint]-b[endpoint]).dropna()
                contrasts[name][endpoint] = dict(n=len(delta),effect=float(delta.mean()) if len(delta) else None)
        joint = all(contrasts[name]['recall_loss']['n']==13 and
                    contrasts[name]['fpr_loss']['n']==23 and
                    contrasts[name]['recall_loss']['effect']>0 and
                    contrasts[name]['fpr_loss']['effect']<=0 for name in BASELINES)
        write_json(dict(primary=METHOD, profile='unchanged hd_posterior', methods=124,
            original_methods_exact=123, calibration_patients=24, evaluation_patients=24,
            threshold=float(threshold), independent_fraction_threshold_exact=True,
            fixed_gate_comparators=BASELINES, point_gate=joint, contrasts=contrasts,
            independent_evaluation=False, production_changed=False),args.out/'report.json')
        for item in rec['input']['files'].values():
            assert file_record(item['path'])['sha256']==item['sha256']
        finish_run(args.out, rec, backend='CPU frozen bulk patient-level loss scorer')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    if jax.default_backend()!='cpu':
        raise RuntimeError('Actual CPU backend required')
    path = Path(__file__).with_name('reference_bagged_profile.py')
    spec = importlib.util.spec_from_file_location('fixed_loss_rank_load', path)
    loader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loader)
    base = loader.load('profile')
    audit = loader.load('bulk_posterior_stage_audit')
    args.out.mkdir(exist_ok=False)
    root = args.task.parent.parent
    rna = argparse.Namespace(out=args.out/'rna', task=args.task)
    scored = argparse.Namespace(out=args.out/'score', task=args.task)
    for part in (rna, scored):
        part.out.mkdir()
    make_rna(rna, root, base)
    score(scored, root, base, audit)


if __name__=='__main__':
    main()
