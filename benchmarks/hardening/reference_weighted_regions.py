"""Research reference-weighted gene RNA effects on fixed native regions."""
import argparse
from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtri
from scipy.stats import pearsonr
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def load(name):
    spec = importlib.util.spec_from_file_location('reference_regions_'+name, Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def weighted_median(values, weights):
    """Exact dyadic weight comparisons, including half-mass ties."""
    values, weights = np.asarray(values, float), np.asarray(weights, float)
    if (values.ndim != 1 or weights.shape != values.shape or not len(values)
            or not np.isfinite(values).all() or not np.isfinite(weights).all() or np.any(weights <= 0)):
        raise ValueError('Finite effects and strictly positive aligned weights required')
    order = np.argsort(values, kind='stable')
    values, weights = values[order], weights[order]
    dyadic = [float(w).as_integer_ratio() for w in weights]
    exponent = max(d.bit_length()-1 for _, d in dyadic)
    integers = [n << (exponent-(d.bit_length()-1)) for n, d in dyadic]
    total, cumulative = sum(integers), 0
    for j, weight in enumerate(integers):
        cumulative += weight
        if 2*cumulative >= total:
            return float((values[j]+values[j+1])/2) if 2*cumulative == total else float(values[j])
    raise ArithmeticError('Weighted median mass was not conserved')


def gene_variance(spread, expected, weighting):
    """Positive descriptive weighting scale, without an inferential SE claim."""
    spread, expected = np.asarray(spread, float), np.asarray(expected, float)
    if (spread.shape != expected.shape or not np.isfinite(spread).all() or np.any(spread < 0)
            or not np.isfinite(expected).all() or np.any(expected < 0)
            or weighting not in ('reference', 'count')):
        raise ValueError('Aligned finite reference spread and supported expected counts required')
    return (spread if weighting == 'reference' else 0)+1/(expected+.5)


def normal_centered_effect(values, normal_effects, gene_bins):
    """Relative RNA against typical normals, with the existing bin-median centre."""
    values, normal_effects, gene_bins = map(np.asarray, (values, normal_effects, gene_bins))
    if (values.ndim != 1 or normal_effects.ndim != 2 or normal_effects.shape[1] != len(values)
            or not len(normal_effects) or gene_bins.shape != values.shape or gene_bins.dtype.kind not in 'iu'
            or not np.array_equal(np.unique(gene_bins), np.arange(int(gene_bins.max())+1))
            or not np.isfinite(values).all() or not np.isfinite(normal_effects).all()):
        raise ValueError('Finite aligned query/normal gene effects and complete integer bin IDs required')
    location = np.median(normal_effects, axis=0)
    effect = values-location
    centre = float(np.median([np.median(effect[gene_bins == j]) for j in np.unique(gene_bins)]))
    return effect-centre, location, centre


def region_locations(values, variance, gene_bins, path, reset):
    values, variance, gene_bins, path, reset = map(np.asarray, (values, variance, gene_bins, path, reset))
    if (values.ndim != 1 or variance.shape != values.shape or gene_bins.shape != values.shape
            or path.ndim != 1 or reset.shape != path.shape or not len(path) or not reset[0]
            or not np.isfinite(values).all() or not np.isfinite(variance).all() or np.any(variance <= 0)
            or not np.array_equal(np.unique(gene_bins), np.arange(len(path)))):
        raise ValueError('Complete gene/bin identities, finite effects and positive variances required')
    group = np.cumsum(np.r_[True, reset[1:] | (path[1:] != path[:-1])])-1
    result = {k:np.empty(len(path)) for k in ('reference_weighted_median', 'reference_weighted_mean', 'unweighted_region_mean')}
    for label in np.unique(group):
        bins = np.flatnonzero(group == label)
        ix = np.flatnonzero(np.isin(gene_bins, bins))
        y, v = values[ix], variance[ix]
        w = 1/v
        if not np.isfinite(w).all():
            raise ArithmeticError('Unresolved numerical precision weights')
        result['reference_weighted_median'][bins] = weighted_median(y, 1/np.sqrt(v))
        result['reference_weighted_mean'][bins] = math.fsum(map(float, w*y))/math.fsum(map(float, w))
        result['unweighted_region_mean'][bins] = math.fsum(map(float, y))/len(y)
    return result


def checks():
    query = np.array([2., 4., 8., 16.]); normals = np.array([[0., 1., 2., 4.], [2., 3., 4., 6.], [4., 5., 6., 8.]])
    bin_ids = np.array([0, 0, 1, 1])
    adjusted, location, centre = normal_centered_effect(query, normals, bin_ids)
    def scalar_median(x):
        ordered = sorted(map(float, x)); n = len(ordered)
        return ordered[n//2] if n%2 else (ordered[n//2-1]+ordered[n//2])/2
    independent_location = [scalar_median(normals[:, j]) for j in range(4)]
    residual = [float(query[j])-independent_location[j] for j in range(4)]
    independent_centre = scalar_median([scalar_median(residual[:2]), scalar_median(residual[2:])])
    np.testing.assert_array_equal(location, independent_location)
    np.testing.assert_array_equal(adjusted, np.array(residual)-independent_centre)
    assert centre == independent_centre
    shift = np.array([8., -4., 16., -8.]); permutation = np.array([2, 0, 3, 1])
    np.testing.assert_array_equal(normal_centered_effect(query+shift, normals+shift, bin_ids)[0], adjusted)
    np.testing.assert_array_equal(normal_centered_effect(query[permutation], normals[:, permutation], bin_ids[permutation])[0], adjusted[permutation])
    np.testing.assert_array_equal(gene_variance([.25, .5], [1.5, 3.5], 'count'), [.5, .25])
    np.testing.assert_array_equal(gene_variance([.25, .5], [1.5, 3.5], 'reference'), [.75, .75])
    # Independent Fraction mass enumeration, with a median interval at half mass.
    fixtures = [([0., 4., 9.], [1., 1., 2.], 6.5), ([0., 4., 9.], [1., 2., 1.], 4.),
                ([0., 4., 9.], [1., np.nextafter(1., 2.), 2.], 4.)]
    for y, w, expected in fixtures:
        fractions = list(map(Fraction, w))
        total = sum(fractions)
        candidates = [j for j in range(len(y)) if sum(fractions[:j+1])*2 >= total]
        j = candidates[0]
        independent = (y[j]+y[j+1])/2 if sum(fractions[:j+1])*2 == total else y[j]
        assert independent == expected == weighted_median(y, w)
        assert weighted_median(y[::-1], w[::-1]) == expected
        assert weighted_median(y, np.array(w)*8) == expected
    values, var, bins = np.array([0., 2., 4., 20., 10., 12.]), np.array([1., 4., 1., 4., 1., 4.]), np.array([0, 0, 0, 1, 2, 3])
    result = region_locations(values, var, bins, np.array([0, 0, 0, 0]), np.array([True, False, True, False]))
    for ix in (np.arange(4), np.arange(4, 6)):
        w = [Fraction(1, int(var[i])) for i in ix]
        exact = float(sum(a*Fraction(float(values[i])) for i, a in zip(ix, w))/sum(w))
        observed = result['reference_weighted_mean'][bins[ix[0]]]
        assert observed == exact
    permutation = np.array([3, 0, 5, 2, 1, 4])
    other = region_locations(values[permutation], var[permutation], bins[permutation], np.array([0, 0, 0, 0]), np.array([True, False, True, False]))
    for name in result:
        np.testing.assert_array_equal(result[name], other[name])
    return dict(exact_fraction_weighted_medians=True, adjacent_float_half_mass_resolved=True,
        exact_fraction_means=True, gene_permutation_and_gap_resets=True, weighting_component_ablation=True,
        independent_normal_gene_location=True, reference_translation_invariance=True)


def panel(args):
    gaussian_regions = getattr(args, 'reference_gaussian_regions', False)
    posterior_boundaries = getattr(args, 'posterior_boundaries', False)
    gaussian_centre = getattr(args, 'gaussian_centre', 'zero')
    gaussian_variance = getattr(args, 'gaussian_variance', 'additive')
    if gaussian_variance not in ('additive','noise-subtracted'):
        raise ValueError('Unsupported Gaussian variance decomposition')
    if gaussian_centre not in ('zero', 'bin-median'):
        raise ValueError('Unsupported Gaussian boundary intercept')
    if gaussian_regions and not args.normal_centering:
        raise ValueError('Reference Gaussian regions require the declared normal-centered gene effects')
    root = args.task.parents[1]
    prior = root/'dev/programme_gene_profile_20261109'/f'refit_{args.reference_refit:03d}_v1'
    base, regions = load('profile'), load('programme_gene_profile')
    manifest = json.loads((prior/'run_manifest.json').read_text())
    native = Path(manifest['input']['files']['genes.tsv.gz']['path']).parent
    source = manifest['input']['files']
    previous = (args.comparison_root/f'refit_{args.reference_refit:03d}_v1'
                if args.comparison_root is not None else prior)
    if args.normal_centering:
        if args.weighting != 'reference' or args.comparison_root is None:
            raise ValueError('Normal centering requires reference weighting and the full24-method comparison root')
        if args.reference_refit == 0:
            previous = args.comparison_root/'bulk_v1'
    if args.weighting == 'count' and args.comparison_root is None:
        raise ValueError('Count-only ablation must retain the full21-method preceding comparator family')
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        scorer=Path(base.__file__), transitive_scorer=Path(base.subarm.__file__),
        region_reference=Path(regions.__file__), prior_manifest=prior/'run_manifest.json',
        previous_predictions=prior/'RNA_predictions.npz', windows=previous/'window_scores.npz',
        previous_scores=previous/'patient_scores.tsv',
        gene_map=root/f'dev/bulk_amplitude_reference_20261109/refit_{args.reference_refit:03d}_v1/gene_map.tsv')
    for name in ('counts', 'samples', 'positions'):
        inputs[name] = Path(source[name]['path'])
    for name in ('genes.tsv.gz', 'bins.tsv', 'samples.tsv'):
        inputs[name] = native/name
    for i in range(68):
        inputs[f'context_{i:03d}'] = prior/f'context_{i:03d}.npz'
    split_path = root/'dev/hardening_20261106_2119/profile_v6/split_manifest.tsv'
    inputs['split'] = split_path
    if gaussian_regions:
        inputs['gaussian_boundary_model'] = Path(__file__).with_name('reference_gaussian_regions.py')
        inputs['native_parameters'] = Path(source['run.json']['path'])
        if posterior_boundaries:
            inputs['posterior_region_average'] = Path(__file__).with_name('posterior_region_average.py')
        if gaussian_centre=='bin-median' or gaussian_variance!='additive' or posterior_boundaries:
            control_root = (root/'dev/reference_gaussian_regions_20261109/bulk_v1'
                if gaussian_variance!='additive' or posterior_boundaries else previous)
            if args.reference_refit and (gaussian_variance!='additive' or posterior_boundaries):
                control_root = root/'dev/reference_gaussian_regions_20261109'/f'refit_{args.reference_refit:03d}_v1'
            inputs['previous_gaussian_decisions'] = control_root/'Gaussian_decisions.json'
            for i in range(68):
                if (control_root/f'Gaussian_context_{i:03d}.npz').exists():
                    inputs[f'previous_gaussian_{i:03d}'] = control_root/f'Gaussian_context_{i:03d}.npz'
        if gaussian_variance!='additive':
            inputs['variance_decomposition'] = Path(__file__).with_name('normal_noise_decomposition.py')
    run = start_run(args.out, inputs, dict(reference_refit=args.reference_refit, weighting=args.weighting,
        normal_centering=args.normal_centering, reference_gaussian_regions=gaussian_regions,
        gaussian_centre=gaussian_centre, gaussian_variance=gaussian_variance,
        posterior_boundaries=posterior_boundaries),
        dict(scope='Tumour RNA effects only; conditional development regression; no calibrated variance claim'), mode='reference_weighted_regions')
    try:
        numerical = checks()
        genes = pd.read_csv(inputs['genes.tsv.gz'], sep='\t')
        bins = pd.read_csv(inputs['bins.tsv'], sep='\t')
        samples = pd.read_csv(inputs['samples'], sep='\t', dtype=str)
        table = pd.read_csv(inputs['counts'], sep='\t', index_col=0)
        G = table.loc[:, samples['sample']].to_numpy(float)[genes.column.to_numpy(int)]
        np.testing.assert_array_equal(table.index.to_numpy()[genes.column.to_numpy(int)], genes.gene)
        corrected, expected, paths = [], [], []
        for i in range(len(samples)):
            with np.load(prior/f'context_{i:03d}.npz') as z:
                corrected.append(z['corrected_gene']); expected.append(z['expected_gene']); paths.append(z['path'])
        corrected, expected, paths = map(np.asarray, (corrected, expected, paths))
        normal = samples.role.eq('normal').to_numpy()
        assert normal.sum() == 20 and len(samples) == 68
        centre = np.median(corrected[normal], axis=0)
        # A weighting heuristic, not an inferential variance or a confidence SE.
        spread = (np.median(abs(corrected[normal]-centre), axis=0)/ndtri(.75))**2
        from haplodose.bins import chromosome_starts, arm_starts
        reset = chromosome_starts(bins) | arm_starts(bins, genome_build='hg38')
        gaussian_records = []
        if gaussian_regions:
            gaussian = load('reference_gaussian_regions')
            numerical['gaussian_boundaries'] = gaussian.checks()
            training = gaussian.train(corrected[normal], expected[normal], genes.bin.to_numpy(int))
            np.savez_compressed(args.out/'Gaussian_reference.npz', normal_bins=training['normal_bins'],
                reference_variance=training['reference_variance'])
            from haplodose.profile import ProfileParams
            gaussian_params = ProfileParams(**json.loads(inputs['native_parameters'].read_text())['params'])
            if posterior_boundaries:
                numerical['posterior_boundaries'] = load('posterior_region_average').checks()
            if gaussian_centre=='bin-median' or gaussian_variance!='additive' or posterior_boundaries:
                old_decisions = {r['sample']:r for r in json.loads(inputs['previous_gaussian_decisions'].read_text())}
                controls = []
                for i in np.flatnonzero(~normal):
                    effect, _, _ = normal_centered_effect(corrected[i], corrected[normal], genes.bin.to_numpy(int))
                    variance = gene_variance(spread, expected[i], args.weighting)
                    fixed = gaussian.fit(effect, expected[i], variance, genes.bin.to_numpy(int),
                        training['reference_variance'], reset, gaussian_params)
                    with np.load(inputs[f'previous_gaussian_{i:03d}']) as old:
                        for key in old.files:np.testing.assert_array_equal(old[key], fixed[key])
                    previous_decision = old_decisions[str(samples.iloc[i]['sample'])]
                    for key in ('beta','tau','grid_scores','near_tied_grid_indices','score_rounding_budget','calibration'):
                        assert previous_decision[key]==fixed[key], (i,key)
                    controls.append(str(samples.iloc[i]['sample']))
                assert len(controls)==48
                write_json(dict(fixed_zero_contexts_exact=controls,decisions_exact=True), args.out/'fixed_zero_control.json')
            boundary_reference_variance=training['reference_variance']
            if gaussian_variance!='additive':
                decomposition=load('normal_noise_decomposition').decompose(
                    corrected[normal],expected[normal],genes.bin.to_numpy(int),training)
                np.savez_compressed(args.out/'normal_noise_decomposition.npz',**decomposition)
                boundary_reference_variance=decomposition['intrinsic_variance']
        names = (('normal_centered_uniform_median', 'normal_centered_reference_median', 'normal_centered_reference_mean')
            if args.normal_centering else (args.weighting+'_weighted_median', args.weighting+'_weighted_mean', 'unweighted_region_mean'))
        if gaussian_regions:
            names = ('gaussian_region_uniform_median', 'gaussian_region_reference_median', 'gaussian_region_reference_mean')
            if gaussian_centre=='bin-median':
                names = ('centered_gaussian_uniform_median', 'centered_gaussian_reference_median', 'centered_gaussian_reference_mean')
            if gaussian_variance!='additive':
                names = ('decomposed_gaussian_uniform_median', 'decomposed_gaussian_reference_median', 'decomposed_gaussian_reference_mean')
            if posterior_boundaries:
                names = ('averaged_gaussian_uniform_median', 'averaged_gaussian_reference_median', 'averaged_gaussian_reference_mean')
        predictions = {name:np.full((len(samples), len(bins)), np.nan) for name in names}
        weights = []; query_centres = np.full(len(samples), np.nan)
        for i in np.flatnonzero(~normal):
            variance = gene_variance(spread, expected[i], args.weighting)
            effect = corrected[i]
            if args.normal_centering:
                effect, location, query_centres[i] = normal_centered_effect(effect, corrected[normal], genes.bin.to_numpy(int))
            selected_path = paths[i]
            if gaussian_regions:
                fitted_boundary = gaussian.fit(effect, expected[i], variance, genes.bin.to_numpy(int),
                    boundary_reference_variance, reset, gaussian_params, centre_mode=gaussian_centre)
                selected_path = fitted_boundary['path']
                array_keys = ('path', 'posterior', 'bin_values', 'sampling_variance', 'total_variance')
                np.savez_compressed(args.out/f'Gaussian_context_{i:03d}.npz', **{k:fitted_boundary[k] for k in array_keys})
                gaussian_records.append(dict(sample=str(samples.iloc[i]['sample']),
                    **{k:v for k,v in fitted_boundary.items() if k not in array_keys}))
            fitted = region_locations(effect, variance, genes.bin.to_numpy(int), selected_path, reset)
            if args.normal_centering:
                predictions[names[0]][i] = regions.region_medians(effect, genes.bin.to_numpy(int), selected_path, reset)
                predictions[names[1]][i] = fitted['reference_weighted_median']
                predictions[names[2]][i] = fitted['reference_weighted_mean']
                if posterior_boundaries:
                    averaged, mass = load('posterior_region_average').average(effect, variance,
                        genes.bin.to_numpy(int), reset, fitted_boundary, gaussian_params)
                    for k, name in enumerate(names):predictions[name][i] = averaged[k]
                    write_json(mass, args.out/f'run_mass_{i:03d}.json')
            else:
                for name, key in zip(names, ('reference_weighted_median', 'reference_weighted_mean', 'unweighted_region_mean')):
                    predictions[name][i] = fitted[key]
            weights.append(variance)
            # Check the unchanged native median from independently retained contexts.
            with np.load(inputs['previous_predictions']) as z:
                np.testing.assert_array_equal(regions.region_medians(corrected[i], genes.bin.to_numpy(int), paths[i], reset),
                    z['programme_gene_native_regions'][i])
        np.savez_compressed(args.out/'RNA_predictions.npz', **predictions)
        if gaussian_regions:
            write_json(gaussian_records, args.out/'Gaussian_decisions.json')
        np.savez_compressed(args.out/'reference_weights.npz', spread=spread, tumour_variance=np.asarray(weights),
            normal_location=centre, query_centres=query_centres,
            sample_ids=samples['sample'].to_numpy(str), normal=normal, genes=genes.gene.to_numpy(str))
        write_json(dict(predictions=file_record(args.out/'RNA_predictions.npz'), weights=file_record(args.out/'reference_weights.npz'),
            numerical=numerical, normal_rows='No candidate estimate; reference roles only',
            normal_ids=samples.loc[normal,'sample'].tolist()), args.out/'RNA_DECISIONS_LOCKED.json')
        # No DNA outcome is read before all RNA effects and weights are frozen.
        with np.load(inputs['windows']) as z:
            truth = {k:z[k] for k in z.files}
        metadata_names = ('samples', 'common', 'loss', 'neutral', 'gain', 'dna', 'split')
        originals = [name for name in truth if name not in metadata_names]
        assert len(originals) == (36 if posterior_boundaries else (33 if gaussian_variance!='additive' else (30 if gaussian_centre=='bin-median' else (27 if gaussian_regions else (24 if args.normal_centering else (21 if args.weighting == 'count' else 18))))))
        windows = {name:truth[name] for name in originals}
        mapping = pd.read_csv(inputs['gene_map'], sep='\t')
        for name, array in predictions.items():
            aligned = base.subarm.align_rows(array, samples['sample'], truth['samples'])
            mapped = base.subarm.gene_values(aligned, genes.set_index('gene').bin, mapping.gene.to_numpy())
            windows[name] = base.subarm.window_median(np.where(mapping.shared.to_numpy(bool)[:,None], mapped, np.nan),
                mapping.window.to_numpy(int), truth['dna'].shape[0])
            assert np.isfinite(windows[name][truth['common']]).all(), name
        np.savez_compressed(args.out/'window_scores.npz', **truth, **{k:v for k,v in windows.items() if k not in truth})
        split = pd.read_csv(split_path, sep='\t', dtype=str).set_index('sample').loc[truth['samples']]
        assert not split.patient.duplicated().any() and np.array_equal(split.split, truth['split'])
        calibration = np.flatnonzero(split.split.eq('calibration'))
        rows = []
        for name, array in windows.items():
            threshold = base.loss_threshold(array, truth['neutral']&truth['common'], calibration)
            losses = base.score_loss(array, truth['loss'], truth['neutral'], truth['common'], threshold, negative=truth['gain']|truth['neutral'])
            for j, sample in enumerate(truth['samples']):
                x,y = array[truth['common'][:,j],j], truth['dna'][truth['common'][:,j],j]
                correlation = float(pearsonr(x,y).statistic) if len(x)>=3 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                if name in names and len(x)>=3 and np.ptp(y)>0 and not np.isfinite(correlation):
                    raise ArithmeticError('Candidate loses an informative patient')
                rows.append(dict(sample=str(sample), split=split.iloc[j].split, method=name, pearson=correlation, **losses[j]))
        scores = pd.DataFrame(rows)
        scores.to_csv(args.out/'patient_scores.tsv', sep='\t', index=False)
        old = pd.read_csv(inputs['previous_scores'], sep='\t', float_precision='round_trip')
        old_check = scores[scores.method.isin(originals)].reset_index(drop=True)
        # CSV decimal round-tripping is checked at the exact frozen float values.
        old_check.to_csv(args.out/'old_method_replay.tsv', sep='\t', index=False)
        pd.testing.assert_frame_equal(old, pd.read_csv(args.out/'old_method_replay.tsv', sep='\t', float_precision='round_trip'), check_exact=True)
        scores.groupby(['split','method'])[['pearson','recall_loss','precision_loss','fpr_loss']].agg(['mean','count']).to_csv(args.out/'summary.tsv', sep='\t')
        paired = []
        selected = scores[scores.split.eq('evaluation')]
        for metric in ('pearson','recall_loss','precision_loss','fpr_loss'):
            matrix = selected.pivot(index='sample', columns='method', values=metric)
            for candidate in (names if args.normal_centering else names[:2]):
                for baseline in windows:
                    if baseline == candidate:
                        continue
                    delta = (matrix[candidate]-matrix[baseline]).dropna()
                    signed = -delta if metric == 'fpr_loss' else delta
                    paired.append(dict(candidate=candidate, baseline=baseline, metric=metric, patients=len(delta),
                        mean=float(delta.mean()), wins=int((signed>0).sum()), losses=int((signed<0).sum()), ties=int((delta==0).sum())))
        continuation_baselines = ('programme_gene_native_regions', 'naive_own_partition')
        if gaussian_regions:
            continuation_baselines += ('normal_centered_reference_median',)
            if gaussian_centre=='bin-median' or gaussian_variance!='additive' or posterior_boundaries:
                continuation_baselines += ('gaussian_region_reference_median',)
            if gaussian_variance!='additive' or posterior_boundaries:
                continuation_baselines += ('centered_gaussian_reference_median',)
            if posterior_boundaries:
                continuation_baselines += ('decomposed_gaussian_reference_median',)
        continuation = [candidate for candidate in (names if args.normal_centering else names[:2]) if all(next(r['mean'] for r in paired
            if r['candidate']==candidate and r['baseline']==baseline and r['metric']=='pearson')>0
            for baseline in continuation_baselines)]
        write_json(dict(methods=list(windows), checks=numerical, paired=paired,
            uncertainty_continuation_candidates=continuation, independent_evaluation=False,
            production_changed=False, default_promoted=False, scientific_criteria_changed=False), args.out/'report.json')
        finish_run(args.out, run, backend='CPU reference-weighted descriptive bulk profiles on frozen native RNA regions')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--reference-refit', type=int, choices=range(9), default=0)
    parser.add_argument('--checks-only', action='store_true')
    parser.add_argument('--weighting', choices=('reference', 'count'), default='reference')
    parser.add_argument('--comparison-root', type=Path)
    parser.add_argument('--normal-centering', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    if args.checks_only:
        write_json(checks(), args.out/'report.json')
    else:
        panel(args)


if __name__ == '__main__':
    main()
