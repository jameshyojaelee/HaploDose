"""Research bin noise from discovery gene composition and reference gene noise.

This is a conditional RNA working model, not a genomic or calibration claim.
"""
import argparse
from contextlib import contextmanager
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys

import jax
import numpy as np
from scipy.stats import nbinom

from haplodose import clones
from haplodose.io import write_json
from haplodose.provenance import file_record, start_run, finish_run


def load(name):
    path = Path(__file__).with_name(name + '.py')
    spec = importlib.util.spec_from_file_location('query_noise_' + name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bin_noise(totals, gbin, gene_size, bins):
    totals, gbin = np.asarray(totals, float), np.asarray(gbin, int)
    sizes = np.broadcast_to(np.asarray(gene_size, float), totals.shape)
    if (totals.ndim != 1 or gbin.shape != totals.shape or not np.isfinite(totals).all()
            or np.any(totals < 0) or np.any(gbin < 0) or np.any(gbin >= bins)
            or not np.isfinite(sizes).all() or np.any(sizes <= 0)):
        raise ValueError('Aligned nonnegative gene totals and positive sizes required')
    pooled = np.bincount(gbin, weights=totals, minlength=bins)
    if np.any(pooled <= 0):
        raise ValueError('Discovery counts do not identify every within-bin composition')
    fractions = totals / pooled[gbin]
    phi = np.bincount(gbin, weights=fractions**2 / sizes, minlength=bins)
    if not np.isfinite(phi).all() or np.any(phi <= 0):
        raise ValueError('Positive bin dispersion not identified')
    return phi, fractions


def numerical_checks():
    totals = np.array([1., 3., 2., 6.])
    gbin = np.array([0, 0, 1, 1])
    sizes = np.array([2., 4., 5., 10.])
    phi, fractions = bin_noise(totals, gbin, sizes, 2)
    exact = [sum(Fraction(int(totals[g]), int(totals[gbin == j].sum()))**2
                 / Fraction(int(sizes[g])) for g in range(4) if gbin[g] == j)
             for j in range(2)]
    operations = 16 * len(totals) + 32
    eps = np.finfo(float).eps
    budget = operations * eps / (1 - operations * eps)
    assert np.max(np.abs(phi - np.array(list(map(float, exact))))) <= budget
    order = np.array([3, 0, 2, 1])
    other, _ = bin_noise(totals[order], gbin[order], sizes[order], 2)
    np.testing.assert_array_equal(phi, other)
    scaled, _ = bin_noise(totals * np.array([2., 2., 4., 4.]), gbin, sizes, 2)
    np.testing.assert_array_equal(phi, scaled)
    # Independent scalar sum of individual gene variances, and SciPy NB moments.
    rows = []
    for j, mean in enumerate((12., 30.)):
        mu = mean * fractions[gbin == j]
        independent = sum(float(m + m*m/r) for m, r in zip(mu, sizes[gbin == j]))
        expected = mean + phi[j]*mean*mean
        scipy_mean, scipy_var = nbinom.stats(1/phi[j], 1/(1+phi[j]*mean), moments='mv')
        allowance = budget * max(1., independent)
        assert abs(independent - expected) <= allowance
        assert abs(float(scipy_mean) - mean) <= allowance
        assert abs(float(scipy_var) - independent) <= allowance
        # Gene-thinned halves have phi/share when gene fractions are held fixed.
        half_var = sum(float(m/2 + m*m/(2*r)) for m, r in zip(mu, sizes[gbin == j]))
        assert abs(half_var - (mean/2 + (phi[j]/.5)*(mean/2)**2)) <= allowance
        rows.append(dict(bin=j, exact_phi=str(exact[j]), phi=float(phi[j]),
                         variance=independent, rounding_allowance=allowance))
    # Confirmation-only perturbation cannot alter this discovery-only function.
    confirmation = np.array([9., 1., 8., 2.])
    before, _ = bin_noise(totals, gbin, sizes, 2)
    confirmation[:] = confirmation[::-1] * 100
    after, _ = bin_noise(totals, gbin, sizes, 2)
    np.testing.assert_array_equal(before, after)
    rejected = False
    try:
        bin_noise(np.array([0., 0., 2., 6.]), gbin, sizes, 2)
    except ValueError:
        rejected = True
    assert rejected
    return dict(passed=True, checks=rows, permutation=True, regional_scale_invariance=True,
                unsupported_bins_rejected=True, confirmation_values_unused=True,
                scope='Algebra and implementation only; not uncertainty or biological validation')


@contextmanager
def membership_noise(registry, records, directory, *, apply):
    original_init = clones.SpecimenModel.__init__
    original_update = clones.SpecimenModel._update_phi
    resolved = {}

    def init(self, *args, **kwargs):
        original_init(self, *args, **kwargs)
        if id(self.rows) not in registry:
            raise RuntimeError('No matched gene fold registry')
        if self.alleles is not None or self.graph is not None or not self.thin:
            raise ValueError('Expression-only, nonspatial gene thinning required')
        info = registry[id(self.rows)]
        # Generate exactly the same cached gene split as the control.
        self.rows.take(np.arange(self.rows.n), fold='prop')
        totals = np.asarray(info['prop'].sum(axis=0)).ravel()
        phi, fractions = bin_noise(totals, info['gbin'], info['resolved_gene_size'], self.J)
        path = directory / f'model_{len(resolved):03d}.npz'
        np.savez_compressed(path, discovery_gene_totals=totals, fractions=fractions,
                            gbin=info['gbin'], gene_size=info['resolved_gene_size'],
                            initial_native_phi=self.phi, proposed_phi=phi, rows=self.rows.rows)
        resolved[id(self)] = dict(phi=phi, index=len(records))
        records.append(dict(kind='initial', apply=apply, arrays=str(path),
            case_seed=info['case_seed'], cells=self.rows.n, bins=self.J,
            native_phi_median=float(np.median(self.phi)), proposed_phi_median=float(np.median(phi)),
            native_over_proposed_median=float(np.median(self.phi/phi)),
            within_bin_composition_source='proposal_gene_counts_only',
            gene_size_source='five_depth_strata_fit_references',
            confirmation_counts_used_for_composition=0))
        if apply:
            self.phi = phi.copy()

    def update(self, R, doses):
        state = resolved[id(self)]
        before = self.phi.copy()
        if apply:
            np.testing.assert_array_equal(self.phi, state['phi'])
        else:
            original_update(self, R, doses)
        records.append(dict(kind='update', apply=apply, model=state['index'],
            before_median=float(np.median(before)), after_median=float(np.median(self.phi)),
            proposed_median=float(np.median(state['phi'])), components=R.shape[1]))

    clones.SpecimenModel.__init__, clones.SpecimenModel._update_phi = init, update
    try:
        yield
    finally:
        clones.SpecimenModel.__init__, clones.SpecimenModel._update_phi = original_init, original_update


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if jax.default_backend() != 'cpu':
        raise RuntimeError('Actual CPU backend required')
    args.out.mkdir(exist_ok=False)
    helper, noise = load('shared_profile_merge_probe'), load('reference_gene_noise')
    modules = ['shared_profile_merge_probe', 'gene_thinning_hook', 'membership_stage_panel', 'membership_weighted_stages',
               'gene_veto_weights', 'recovery_baseline', 'reference_gene_noise']
    run = start_run(args.out, dict(script=Path(__file__), design=args.task/'DESIGN.md',
        **{name: Path(__file__).with_name(name+'.py') for name in modules}),
        dict(arms=['baseline', 'sham', 'candidate'], depth_strata=5),
        dict(scope='Exposed RNA recovery mechanism, biological n=0; no public promotion'),
        mode='query_composition_membership_noise')
    original_load, argv = noise.load, sys.argv
    traces, summaries = {}, {}
    try:
        write_json(numerical_checks(), args.out/'numerics.json')
        for arm in ('baseline', 'sham', 'candidate'):
            directory = args.out/(arm+'_noise')
            directory.mkdir()
            records = []

            def wrapped_load(name):
                module = original_load(name)
                if name == 'gene_thinning_hook' and arm != 'baseline':
                    original = module.gene_native

                    @contextmanager
                    def wrapped(*a, **k):
                        with original(*a, **k) as registry:
                            with membership_noise(registry, records, directory, apply=arm == 'candidate'):
                                yield registry
                    module.gene_native = wrapped
                return module

            noise.load = wrapped_load
            sys.argv = [str(noise.__file__), '--out', str(args.out/arm), '--task', str(args.task),
                        '--stage', 'both', '--depth-strata', '5']
            noise.main()
            summaries[arm] = json.loads((args.out/arm/'report.json').read_text())
            traces[arm] = records
            write_json(records, args.out/(arm+'_noise_trace.json'))
            if arm == 'sham':
                exact = helper.exact_outputs(args.out/'baseline/panel', args.out/'sham/panel')
                assert len(exact) == 5 and records
                write_json(dict(passed=True, cases=exact), args.out/'SHAM_LOCKED.json')
            if arm != 'baseline':
                for path in (args.out/'baseline/folds').glob('*.npz'):
                    with np.load(path) as old, np.load(args.out/arm/'folds'/path.name) as new:
                        keys = [key for key in old.files if not key.startswith('veto')]
                        for key in keys:
                            assert key in new.files
                            np.testing.assert_array_equal(old[key], new[key])
        before, after = summaries['baseline']['cases'], summaries['candidate']['cases']
        comparisons = []
        for old, new in zip(before, after):
            assert old['seed'] == new['seed']
            comparisons.append(dict(seed=old['seed'], before=old, after=new,
                criterion_losses=[k for k, v in old['criteria'].items()
                                  if v['passed'] and not new['criteria'][k]['passed']],
                added_s3_groups=new['s3_groups'] > old['s3_groups']))
        advance = (len(comparisons) == 5 and all(r['after']['original_assertion_passed'] for r in comparisons)
                   and not any(r['criterion_losses'] or r['added_s3_groups'] for r in comparisons))
        write_json(dict(baseline_recovered=summaries['baseline']['recovered'],
            candidate_recovered=summaries['candidate']['recovered'], comparisons=comparisons,
            continuation_passed=advance, sham_exact=True, gene_folds_exact=True,
            original_criteria_unchanged=True, production_changed=False, biological_specimens=0,
            false_split_control_established=False, promotion=False), args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU installed wheel, research-only composition noise')
    except BaseException as error:
        for arm, records in traces.items():
            write_json(records, args.out/(arm+'_noise_trace.json'))
        finish_run(args.out, run, error=error)
        raise
    finally:
        noise.load, sys.argv = original_load, argv


if __name__ == '__main__':
    main()
