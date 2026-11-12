"""Exact rational reference for the frozen 27 paired discordance comparisons."""
import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path

from haplodose.provenance import start_run, finish_run, file_record


def exact_tail(lost, gained):
    n = lost + gained
    return min(Fraction(1), Fraction(2 * sum(math.comb(n, j)
        for j in range(min(lost, gained) + 1)), 2 ** n))


def exact_holm(values):
    # Independent unsorted formula, including tied raw probabilities.
    return [min(Fraction(1), max((len(values) - sum(v < u for v in values)) * u
        for u in values if u <= p)) for p in values]


def reference_checks():
    # Enumerate the null's equiprobable binary assignments, without tail sums.
    for n in range(11):
        totals = [sum(bits) for bits in itertools.product((0, 1), repeat=n)]
        for gained in range(n + 1):
            extreme = sum(abs(2 * v - n) >= abs(2 * gained - n) for v in totals)
            if exact_tail(n - gained, gained) != Fraction(extreme, 2 ** n):
                raise AssertionError('Exact tail disagrees with exhaustive assignments')
    example = [Fraction(1, 100), Fraction(1, 100), Fraction(1, 5), Fraction(1)]
    if exact_holm(example) != [Fraction(1, 25), Fraction(1, 25), Fraction(2, 5), Fraction(1)]:
        raise AssertionError('Tied Holm fixture failed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    inputs = dict(script=Path(__file__), summary=args.task/'summary_v1/summary.json',
        decisions=args.task/'summary_v1/decisions.json', audit=args.task/'audit_v1/report.json',
        audit_manifest=args.task/'audit_v1/run_manifest.json')
    manifest = json.loads(inputs['audit_manifest'].read_text())
    if manifest['status'] != 'complete':
        raise AssertionError('Primary independent audit incomplete')
    expected = manifest['output']['files']['report.json']
    if file_record(inputs['audit'])['sha256'] != expected['sha256']:
        raise AssertionError('Primary audit changed')
    args.out.mkdir(exist_ok=False)
    record = start_run(args.out, inputs, dict(family=27, alpha='1/20'),
        dict(scope='Numerical reference only; no new scientific endpoint or selection'),
        mode='joint_paired_reference')
    try:
        reference_checks()
        summary = json.loads(inputs['summary'].read_text())
        decisions = json.loads(inputs['decisions'].read_text())
        contrasts = summary['paired']
        if len(contrasts) != 27:
            raise AssertionError('Paired family changed')
        exact = []
        for contrast in contrasts:
            values = []
            for arm in ('joint_proposal', 'joint_full'):
                selected = [r for r in decisions if r['arm'] == arm and
                    all(r[key] == contrast[key] for key in ('condition', 'method', 'alpha'))]
                if sorted(r['replicate'] for r in selected) != list(range(200)):
                    raise AssertionError('Paired population changed')
                values.append({r['replicate']: r['event'] for r in selected})
            lost = sum(values[0][i] and not values[1][i] for i in range(200))
            gained = sum(values[1][i] and not values[0][i] for i in range(200))
            if (lost, gained) != (contrast['lost'], contrast['gained']):
                raise AssertionError('Discordance counts disagree')
            exact.append(exact_tail(lost, gained))
        adjusted = exact_holm(exact)
        rows = []
        for row, raw, holm in zip(contrasts, exact, adjusted):
            # The exact rational calculation determines the scientific decision.
            # Retain all binary floating-point discrepancies without rounding them
            # away or substituting a permissive agreement tolerance.
            candidate = dict(row, exact_nominal=str(raw), exact_holm=str(holm),
                nominal_reference=float(raw), holm_reference=float(holm),
                nominal_discrepancy=row['nominal_p'] - float(raw),
                holm_discrepancy=row['holm_p'] - float(holm),
                exact_holm_reject=holm <= Fraction(1, 20),
                recorded_holm_reject=row['holm_p'] <= .05)
            rows.append(candidate)
        (args.out/'report.json').write_text(json.dumps(dict(contrasts=rows,
            complete_paired_populations=27, exhaustive_null_max_n=10,
            scientific_thresholds_changed=False), indent=2))
        if any(r['exact_holm_reject'] != r['recorded_holm_reject'] for r in rows):
            raise AssertionError('Exact rational Holm decision differs')
        finish_run(args.out, record, backend='cpu')
    except Exception as error:
        finish_run(args.out, record, error=error)
        raise


if __name__ == '__main__':
    main()
