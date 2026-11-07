"""Small observational-equivalence and input-contract falsifiers; compute node only.

This diagnostic uses synthetic measurements, not biological replication. Its
exact equivalences delimit possible claims; they do not rank predictive methods.
All files are created below the explicitly supplied task root and fresh output.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, special, stats

from haplodose.bulk import run_bulk
from haplodose.clones import mixture_weights
from haplodose.io import new_output_dir, write_json
from haplodose.kernels.estep import cell_loglik
from haplodose.profile import ProfileParams, fit_pooled, segment
from haplodose.provenance import file_record, finish_run, start_run


def fresh_output(task_root, output):
    root, output = task_root.resolve(strict=True), output.resolve()
    relative = output.relative_to(root)
    if not relative.parts or relative.parts[0] == 'baseline':
        raise ValueError('Probe output must be a new descendant of task_root outside baseline')
    if output.exists():
        raise FileExistsError(output)
    return new_output_dir(output)


def composition(values):
    return np.asarray(values, float) / np.sum(values)


def likelihood_contrast(y, p, q, phi):
    """Independent float64 NB log-likelihood contrast, conditional working model."""
    L = y.sum()
    first = stats.nbinom.logpmf(y, 1 / phi, 1 / (1 + phi * L * p))
    second = stats.nbinom.logpmf(y, 1 / phi, 1 / (1 + phi * L * q))
    return float(np.sum(second - first))


def equivalents(rng, J):
    base = rng.dirichlet(np.ones(J) * 3)
    region = np.arange(J) < J // 4
    copy3, copy4 = np.ones(J) * 2, np.ones(J) * 2
    copy3[region], copy4[region] = 3, 4
    first_purity = composition(base * (1 - .5 + .5 * copy3 / 2))
    second_purity = composition(base * (1 - .25 + .25 * copy4 / 2))
    dosage = np.where(region, .6 * np.log(1.5), 0.)
    rate = composition(base * np.exp(dosage))
    examples = [
        ('global_log_rate_shift', rate, special.softmax(np.log(rate) + .7),
         'global RNA amount is unidentified after compositional normalization'),
        ('purity_copy_number', first_purity, second_purity,
         'purity=.5,CN=3 versus purity=.25,CN=4 under linear dosage response'),
        ('balanced_genome_doubling', base, composition(2 * base),
         'balanced 1:1 versus 2:2 alleles also both have BAF=.5'),
        ('transcription_versus_dosage', rate, composition(base * np.exp(dosage)),
         'regional transcription without a DNA event can reproduce the RNA dosage mean')]
    phi = np.linspace(.05, .5, J)
    rows = []
    for name, p, q, explanation in examples:
        y = rng.multinomial(20000, p)
        package = cell_loglik(y[None], np.array([y.sum()]), np.stack([np.log(p), np.log(q)]), phi)[0]
        # Standard floating-point reduction budgets, not scientific tolerances.
        eps = np.finfo(np.float32).eps
        budget = (J + 16) * eps / (1 - (J + 16) * eps) * (y.sum() + np.sum(1 / phi))
        rows.append(dict(case=name, max_rate_difference=float(np.max(np.abs(p - q))),
                         scipy_float64_nb_contrast=likelihood_contrast(y, p, q, phi),
                         package_float32_nb_contrast=float(package[1] - package[0]),
                         arithmetic_budget=float(budget), interpretation=explanation,
                         estimand='conditional_relative_RNA_composition'))
    return rows


def reference_absorption(rng, J):
    base = np.full(J, 1 / J)
    altered = composition(base * np.exp(np.r_[np.full(J // 4, .6 * np.log(1.5)), np.zeros(3 * J // 4)]))
    references = np.stack([base, altered], axis=1)
    y = rng.multinomial(500000, altered).astype(float)
    fitted_weights, mixture_diagnostic = mixture_weights(y, references, return_diagnostics=True)
    objective = lambda w: -float(y @ np.log(references @ np.array([1 - w, w])))
    independent = optimize.minimize_scalar(objective, bounds=(0, 1), method='bounded',
                                          options={'xatol': np.sqrt(np.finfo(float).eps)})
    candidates = np.array([0., independent.x, 1.])
    optimum = float(candidates[np.argmin([objective(w) for w in candidates])])
    independent_gap = float(objective(fitted_weights[1]) - objective(optimum))
    objective_arithmetic = (J + 16) * np.finfo(float).eps * abs(objective(optimum))
    gap_budget = y.sum() * mixture_diagnostic['float64_roundoff_budget'] + objective_arithmetic
    sanity_gate = bool(mixture_diagnostic['converged'] and
        mixture_diagnostic['normalized_objective_gap_upper'] <= mixture_diagnostic['float64_roundoff_budget'] and
        abs(independent_gap) <= gap_budget)
    if not sanity_gate:
        raise RuntimeError(f'Reference mixture failed independent Poisson-MLE sanity gate: gap={independent_gap}, budget={gap_budget}')
    reset = np.zeros(J, bool)
    reset[::J // 4] = True
    params = ProfileParams(beta=.6, beta_bounds=(.6, .6))
    fits = []
    for name, weights in [('fixed_correct_identity', np.array([1., 0.])),
                          ('freely_fitted_reference', fitted_weights),
                          ('fixed_alternative_identity', np.array([0., 1.]))]:
        E = y.sum() * (references @ weights)
        center, beta, tau = fit_pooled(y, E, reset, params)
        post, path, dose, _ = segment(y[None], E[None], reset, [center], beta, [tau], params)
        fits.append(dict(setting=name, reference_weights=weights.tolist(), center=center, beta=beta, tau=tau,
                         regional_mean_rna_effect=float(dose[0, :J // 4].mean()),
                         regional_gain_fraction=float((path[0, :J // 4] == 2).mean()),
                         regional_neutral_fraction=float((path[0, :J // 4] == 1).mean())))
    return dict(fits=fits, independent_alternative_weight=optimum,
                fitted_objective_gap=independent_gap, objective_gap_arithmetic_budget=gap_budget,
                mixture_optimality=mixture_diagnostic, independent_mle_sanity_gate=sanity_gate,
                profile_params=asdict(params),
                interpretation='Reference composition and a regional dosage effect are confounded in this construction; '
                               'reference identity must be supplied independently to distinguish them')


def count_contracts(out):
    inputs = new_output_dir(out / 'contract_inputs')
    genes = [f'G{i}' for i in range(8)]
    samples = pd.DataFrame({'sample': ['N0', 'N1', 'N2', 'T0'], 'role': ['normal'] * 3 + ['tumour']})
    samples.to_csv(inputs / 'samples.tsv', sep='\t', index=False)
    pd.DataFrame({'gene': genes, 'chromosome': 'chr1', 'start': np.arange(8) * 1000,
                  'end': np.arange(8) * 1000 + 500}).to_csv(inputs / 'positions.tsv', sep='\t', index=False)
    cases = []
    for case, kind in [('undeclared', None), ('declared_TPM', 'TPM'), ('fractional_raw', 'raw'),
                       ('negative_raw', 'raw'), ('zero_library_raw', 'raw')]:
        values = np.full((8, 4), 20.)
        if case == 'fractional_raw':
            values[0, 3] = 20.5
        elif case == 'negative_raw':
            values[0, 3] = -1
        elif case == 'zero_library_raw':
            values[:, 3] = 0
        counts = inputs / (case + '.tsv')
        pd.DataFrame(values, index=genes, columns=samples['sample']).to_csv(counts, sep='\t')
        try:
            run_bulk(counts, inputs / 'samples.tsv', inputs / 'positions.tsv', out / ('reject_' + case),
                     counts_kind=kind, genome_build='hg38', coordinate_system='zero-based-half-open', gene_namespace='custom')
        except ValueError as error:
            cases.append(dict(case=case, rejected=True, reason=str(error)))
        else:
            cases.append(dict(case=case, rejected=False, reason='unexpected acceptance'))
        manifest = out / ('reject_' + case) / 'run_manifest.json'
        if manifest.exists():
            child = json.loads(manifest.read_text())
            cases[-1]['run_manifest'] = file_record(manifest)
            cases[-1]['output_id'] = child.get('output', {}).get('id')
    return dict(cases=cases, all_expected_rejections=all(c['rejected'] for c in cases),
                fixture_files={p.name: file_record(p) for p in sorted(inputs.iterdir())},
                limitation='Integer-valued normalized expression can mimic raw counts numerically; an honest '
                           'measurement declaration and source provenance remain identifying assumptions')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--seed', type=int, default=20261106)
    parser.add_argument('--bins', type=int, default=120)
    args = parser.parse_args()
    if args.bins < 12 or args.bins % 4:
        parser.error('bins must be >=12 and divisible by four')
    out = fresh_output(args.task_root, args.out)
    record = start_run(out, {'probe_source': __file__}, {'seed': args.seed, 'bins': args.bins},
                       {'measurement': 'synthetic counts and exact mean compositions',
                        'replication': 'no biological replication'}, mode='identifiability_diagnostic')
    try:
        rng = np.random.default_rng(args.seed)
        result = dict(equivalences=equivalents(rng, args.bins),
                      reference_absorption=reference_absorption(rng, args.bins),
                      input_contract=count_contracts(out),
                      evidence_scope='computational checks and observational nonidentifiability, not DNA validation')
        write_json(result, out / 'summary.json')
        if not result['input_contract']['all_expected_rejections']:
            raise AssertionError('A malformed measurement contract was accepted')
        if any(abs(row['package_float32_nb_contrast']) > row['arithmetic_budget']
               for row in result['equivalences']):
            raise AssertionError('Equivalent observations disagree beyond the numerical error budget')
        finish_run(out, record, backend='cpu')
        print(json.dumps({'summary': str(out / 'summary.json'), 'count_rejections': result['input_contract']['all_expected_rejections']}))
    except Exception as error:
        finish_run(out, record, error=error)
        raise


if __name__ == '__main__':
    main()
