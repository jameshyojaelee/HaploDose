"""Research transfer of the validated RNA smoother into assignment refinement."""
import argparse
from contextlib import contextmanager
from fractions import Fraction
import hashlib
import importlib.util
import inspect
import itertools
import json
import math
from pathlib import Path
import sys

import jax
import numpy as np

from haplodose import clones, pipeline
from haplodose.io import write_json
from haplodose.provenance import file_record, start_run, finish_run
from haplodose.stacking import stacked_rna_profile


def load(name):
    path = Path(__file__).with_name(name+'.py')
    spec = importlib.util.spec_from_file_location('gene_profile_member_'+name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ordered_profile(previous, proposed, weights, common):
    """Weighted isotonic projection with exact previous equality constraints.

    Preserve old pairwise order and ties; merge adjacent order violations.
    The exposure-weighted mean is the proposed common tumour profile.
    """
    previous, proposed, weights = map(lambda a: np.asarray(a, float), (previous, proposed, weights))
    if (previous.shape != proposed.shape or previous.ndim != 2 or weights.shape != (len(previous),)
            or np.asarray(common).shape != (previous.shape[1],)
            or not all(np.isfinite(a).all() for a in (previous, proposed, weights, common))
            or np.any(weights <= 0)):
        raise ValueError('Aligned finite tumour profiles and positive exposures required')
    result = np.empty_like(previous)
    for j in range(previous.shape[1]):
        groups = [np.flatnonzero(previous[:, j] == value) for value in np.unique(previous[:, j])]
        blocks = []
        for group in groups:
            w = math.fsum(float(weights[k]) for k in group)
            total = math.fsum(float(weights[k]*proposed[k, j]) for k in group)
            blocks.append((group.tolist(), w, total))
            while len(blocks) > 1 and blocks[-2][2]/blocks[-2][1] > blocks[-1][2]/blocks[-1][1]:
                a, b = blocks[-2:]
                blocks[-2:] = [(a[0]+b[0], math.fsum((a[1], b[1])), math.fsum((a[2], b[2])))]
        for group, w, total in blocks:
            result[group, j] = total/w
        centre = math.fsum(float(weights[k]*result[k, j]) for k in range(len(weights)))/math.fsum(map(float, weights))
        result[:, j] = result[:, j]-centre+common[j]
    return result


def checks():
    old = np.array([[0., 0.], [0., 1.], [1., 2.], [2., 3.]])
    raw = np.array([[0., 4.], [2., 3.], [1., -1.], [4., 2.]])
    weight = np.array([1., 1., 2., 1.])
    common = np.array([.25, -.5])
    actual = ordered_profile(old, raw, weight, common)
    references = []
    # Independent exhaustive contiguous partitions of the ordered old levels.
    # Fraction SSE chooses the constrained optimum, rather than repeating PAVA.
    for j in range(2):
        levels = sorted(set(old[:, j]))
        candidates = []
        for mask in range(1 << (len(levels)-1)):
            cuts = [0]+[i+1 for i in range(len(levels)-1) if mask & (1 << i)]+[len(levels)]
            fitted = [None]*len(weight)
            means = []
            for lo, hi in zip(cuts, cuts[1:]):
                indices = [k for k in range(len(weight)) if old[k, j] in levels[lo:hi]]
                denominator = sum(Fraction(int(weight[k])) for k in indices)
                mean = sum(Fraction(int(weight[k]))*Fraction(int(raw[k, j])) for k in indices)/denominator
                means.append(mean)
                for k in indices:
                    fitted[k] = mean
            if any(a > b for a, b in zip(means, means[1:])):
                continue
            loss = sum(Fraction(int(weight[k]))*(fitted[k]-Fraction(int(raw[k, j])))**2 for k in range(len(weight)))
            candidates.append((loss, fitted))
        _, reference = min(candidates, key=lambda a: a[0])
        mean = sum(Fraction(int(weight[k]))*reference[k] for k in range(len(weight)))/sum(map(Fraction, map(int, weight)))
        reference = [v-mean+Fraction(float(common[j])) for v in reference]
        references.append(reference)
        operations = 128*len(weight)+64
        gamma = operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)
        allowance = gamma*max(1., float(np.max(abs(raw))))
        assert np.max(abs(actual[:, j]-np.asarray(list(map(float, reference))))) <= allowance
        for a, b in itertools.product(range(len(weight)), repeat=2):
            if old[a, j] == old[b, j]:
                assert actual[a, j] == actual[b, j]
            if old[a, j] < old[b, j]:
                assert actual[a, j] <= actual[b, j]
    permutation = np.array([3, 0, 2, 1])
    permuted = ordered_profile(old[permutation], raw[permutation], weight[permutation], common)
    np.testing.assert_array_equal(permuted[np.argsort(permutation)], actual)
    np.testing.assert_array_equal(ordered_profile(old, raw, weight*2, common), actual)
    assert np.isfinite(actual).all()
    return dict(passed=True, independent_fraction_partition_optima=[[str(v) for v in x] for x in references],
                exact_ties_preserved=True, pairwise_order_preserved=True, exposure_scale_invariant=True,
                cell_component_permutation=True, scope='Projection arithmetic only; no genomic validation')


@contextmanager
def gene_profiles(records, directory, *, apply):
    original_fit = pipeline.fit_specimen
    original_shared = clones.SpecimenModel.shared_profiles
    signature = inspect.signature(original_fit)
    registry, cache = {}, {}

    def fit(*args, **kwargs):
        bound = signature.bind(*args, **kwargs)
        bound.apply_defaults()
        arguments = bound.arguments
        rows, loader = arguments['rows'], arguments['gene_loader']
        if arguments['alleles'] is not None or arguments['graph'] is not None or loader is None:
            raise ValueError('This transfer experiment requires expression-only nonspatial gene counts')
        frame = sys._getframe(1)
        if frame.f_code is not pipeline._run_single.__code__:
            raise RuntimeError('Missing exact pipeline reference context')
        context = frame.f_locals
        keys = context['type_keys']
        if any(key not in context['ref_gene'] for key in keys):
            raise ValueError('No in-sample or unrecorded reference composition allowed')
        X, gbin = loader(np.arange(rows.n))
        reference = np.vstack([context['ref_gene'][key] for key in keys])
        np.testing.assert_array_equal(gbin, context['used'].bin)
        registry[id(rows)] = dict(rows=rows, X=X.astype(float), gbin=np.asarray(gbin, int),
            references=reference, source=str(context['counts']), genes=context['used'].gene.tolist(),
            cell_ids=context['meta'].cell_id.to_numpy()[rows.rows].tolist(), keys=keys)
        path = directory/f'input_{len(registry)-1:03d}.npz'
        np.savez_compressed(path, reference_gene_counts=reference, gbin=gbin, rows=rows.rows,
                            genes=np.asarray(registry[id(rows)]['genes']),
                            cell_ids=np.asarray(registry[id(rows)]['cell_ids']))
        registry[id(rows)]['record'] = str(path)
        return original_fit(*bound.args, **bound.kwargs)

    def smooth(effects, gbin, blocks):
        key = hashlib.sha256(effects.tobytes()+gbin.tobytes()+blocks.tobytes()).hexdigest()
        if key not in cache:
            profile, diagnostics = stacked_rna_profile(effects, gbin, blocks,
                gene_segments=True, average_segments=True)
            cache[key] = profile, diagnostics
        return cache[key]

    def shared(self, Y, L, XL, labels=None):
        result = original_shared(self, Y, L, XL, labels=labels)
        caller = sys._getframe(1)
        if caller.f_code is not clones.SpecimenModel.em.__code__:
            raise RuntimeError('Assignment replacement outside the declared refinement stage')
        info = registry[id(self.rows)]
        R = caller.f_locals['R']
        genes = np.asarray(R.astype(float).T@info['X'])
        gbin = info['gbin']
        # Verify the actual soft-gene pooling conserves the native bin pools.
        aggregated = np.stack([genes[:, gbin == j].sum(axis=1) for j in range(self.J)], axis=1)
        operations = 32*self.rows.n*len(gbin)+256
        unit = operations*np.finfo(float).eps
        if unit >= 1:
            raise ArithmeticError('Unresolved gene pooling rounding budget')
        budget = unit/(1-unit)*max(1., float(np.max(abs(Y))))
        if np.max(abs(aggregated-Y)) > budget:
            raise ArithmeticError('Gene and native pooled bin counts disagree')
        tumour = np.arange(1, len(Y))
        old = np.asarray(result['dose']).copy()
        if not len(tumour) or np.any(np.asarray(L)[tumour] <= 0):
            records.append(dict(kind='abstention', reason='no_positive_tumour_pool', apply=apply))
            return result
        share = pipeline.gene_shares(info['references'], self.weights, gbin, self.J)
        expected = np.asarray(result['E'])
        effects = np.log((genes+.5)/(expected[:, gbin]*share[None]+.5))
        effects -= np.median(effects, axis=1, keepdims=True)
        blocks = np.cumsum(self.reset).astype(np.int64)
        raw, diagnostics = [], []
        for k in tumour:
            profile, decision = smooth(effects[k], gbin, blocks)
            raw.append(profile)
            diagnostics.append(decision)
        # The common baseline uses exactly the same native pooled expectation
        # and existing reference-gene allocation; no reference is DNA-selected.
        pooled = self.profiles(np.asarray(Y)[tumour].sum(axis=0), np.asarray(L)[tumour].sum(),
                               np.asarray(XL)[tumour].sum(axis=0))
        gene_total = genes[tumour].sum(axis=0)
        common_effect = np.log((gene_total+.5)/(pooled['E'][0][gbin]*share+.5))
        common_effect -= np.median(common_effect)
        common, common_decision = smooth(common_effect, gbin, blocks)
        proposed = old.copy()
        proposed[tumour] = ordered_profile(old[tumour], np.asarray(raw), np.asarray(L)[tumour], common)
        for a, b in itertools.product(tumour, repeat=2):
            tied = old[a] == old[b]
            assert np.array_equal(proposed[a, tied], proposed[b, tied])
            if np.any(proposed[a, old[a] < old[b]] > proposed[b, old[a] < old[b]]):
                raise ArithmeticError('Replacement reverses the existing pairwise order')
        number = len(records)
        path = directory/f'profile_{number:04d}.npz'
        np.savez_compressed(path, R=R, gene_pools=genes, bin_pools=Y, E=expected, L=L, XL=XL,
                            share=share, gbin=gbin, effects=effects, common_effect=common_effect,
                            native_dose=old, proposed_dose=proposed, raw_gene_profiles=np.asarray(raw),
                            common=common, blocks=blocks)
        detail = directory/f'profile_{number:04d}.json'
        write_json(dict(components=diagnostics, common=common_decision), detail)
        records.append(dict(kind='profile', apply=apply, arrays=str(path), diagnostics=str(detail),
            source_counts=info['source'], input=info['record'], cells=self.rows.n,
            components=len(Y), gene_pool_discrepancy=float(np.max(abs(aggregated-Y))),
            pooling_rounding_budget=budget, maximum_dose_change=float(np.max(abs(proposed-old))),
            pairwise_ties_and_order_preserved=True))
        write_json(records, directory/'records.json')
        if apply:
            if 'dose_own' not in result:
                result['dose_own'] = old
            result['dose'] = proposed
        return result

    pipeline.fit_specimen, clones.SpecimenModel.shared_profiles = fit, shared
    try:
        yield
    finally:
        pipeline.fit_specimen, clones.SpecimenModel.shared_profiles = original_fit, original_shared


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if jax.default_backend() != 'cpu':
        raise RuntimeError('Actual CPU backend required')
    args.out.mkdir(exist_ok=False)
    runner, helper = load('membership_stage_panel'), load('shared_profile_merge_probe')
    root = args.task.parents[1]
    panel = root/'dev/recovery_20261107/recovery_panel_v1'
    assertions = root/'dev/recovery_20261107/baseline/tests/test_split_mode.py'
    run = start_run(args.out, dict(script=Path(__file__), design=args.task/'DESIGN.md',
        runner=Path(runner.__file__), helper=Path(helper.__file__), panel=panel/'manifest.json',
        assertions=assertions), dict(arms=['baseline', 'sham', 'candidate'], original_criteria=30),
        dict(scope='Gene RNA profile transfer to refinement; exposed simulations; biological n=0'),
        mode='gene_profile_membership')
    reports, argv = {}, sys.argv
    try:
        write_json(checks(), args.out/'numerics.json')
        for arm in ('baseline', 'sham', 'candidate'):
            directory = args.out/(arm+'_profiles')
            directory.mkdir()
            trace = []
            sys.argv = [str(runner.__file__), '--out', str(args.out/arm), '--panel-root', str(panel),
                        '--assertion-source', str(assertions), '--stage', 'hard']
            if arm == 'baseline':
                runner.main()
            else:
                with gene_profiles(trace, directory, apply=arm == 'candidate'):
                    runner.main()
                assert trace and any(row['kind'] == 'profile' for row in trace)
            reports[arm] = json.loads((args.out/arm/'summary.json').read_text())
            write_json(trace, args.out/(arm+'_profile_trace.json'))
            if arm == 'sham':
                exact = helper.exact_outputs(args.out/'baseline', args.out/'sham')
                assert len(exact) == 5
                write_json(dict(passed=True, cases=exact), args.out/'SHAM_LOCKED.json')
        comparisons = []
        for before, after in zip(reports['baseline']['cases'], reports['candidate']['cases']):
            assert before['seed'] == after['seed']
            comparisons.append(dict(seed=before['seed'], before=before, after=after,
                criterion_losses=[key for key, value in before['criteria'].items()
                                  if value['passed'] and not after['criteria'][key]['passed']],
                additional_s3_groups=after['s3_groups'] > before['s3_groups']))
        advance = (len(comparisons) == 5 and all(r['after']['original_assertion_passed'] for r in comparisons)
                   and not any(r['criterion_losses'] or r['additional_s3_groups'] for r in comparisons))
        write_json(dict(comparisons=comparisons,
            baseline_recovered=sum(r['original_assertion_passed'] for r in reports['baseline']['cases']),
            candidate_recovered=sum(r['original_assertion_passed'] for r in reports['candidate']['cases']),
            continuation_passed=advance, original_criteria_unchanged=True, biological_specimens=0,
            production_changed=False, false_split_control_established=False, promotion=False), args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU installed profile wheel; research assignment smoother')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise
    finally:
        sys.argv = argv


if __name__ == '__main__':
    main()
