"""Independent counterexample and old/new parity for exact confirmation ties."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy import stats

from haplodose.confirmation import Contrast, FrozenDiscovery, score_confirmation
from haplodose.provenance import start_run, finish_run, file_record


def fixtures():
    counts = np.array([[a, b, 13, 7] for a in range(1, 11) for b in range(1, 11)])
    for name, rates, prior, weights in (
        ('symmetric', [[.4, .1, .4, .1], [.1, .4, .1, .4]], [.5, .5], [[4, 1], [1, 4]]),
        ('identical_proposal', [[.4, .1, .4, .1], [.4, .1, .1, .4]], [.5, .5], [[4, 1], [4, 1]]),
        ('unequal_prior', [[.4, .1, .4, .1], [.4, .1, .1, .4]], [.25, .75], [[4, 1], [4, 1]]),
    ):
        artifact = FrozenDiscovery(log_rates=np.log(rates), phi=np.full(4, .1),
            proposal_mask=[True, True, False, False], test_mask=[False, False, True, True],
            contrasts=(Contrast(0, 1, (2,), 1, 'region'),), discovery_cell_ids=('train',),
            candidate_prior=prior)
        yield name, counts, artifact, weights


def exact_candidates(counts, weights, priors):
    """Rational NB probabilities: integer size 10, rational rates and priors."""
    labels = []
    for row in counts:
        probabilities = []
        depth = int(sum(row[:2]))
        for weight, prior in zip(weights, priors):
            value = Fraction(float(prior))
            for count, numerator in zip(row[:2], weight):
                count = int(count)
                mean = Fraction(depth * numerator, sum(weight))
                p = 10 / (10 + mean)
                value *= math.comb(count + 9, count) * p**10 * (1-p)**count
            probabilities.append(value)
        winners = [k for k, value in enumerate(probabilities) if value == max(probabilities)]
        labels.append(winners[0] if len(winners) == 1 else -1)
    return np.array(labels)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--reference', type=Path)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
                  baseline=args.task/'baseline/source_manifest.json')
    if args.reference:
        inputs['old_results'] = args.reference/'report.json'
        inputs['old_manifest'] = args.reference/'run_manifest.json'
    run = start_run(args.out, inputs, dict(exact_ties=True),
                    dict(scope='Computational identity, not biological calibration'), mode='confirmation_ties')
    try:
        old = json.loads(inputs['old_results'].read_text()) if args.reference else None
        records = {}
        for name, Y, artifact, weights in fixtures():
            ids = [f'cell{i}' for i in range(len(Y))]
            result = score_confirmation(Y, ids, artifact)
            scores = result['proposal_log_scores']
            # Independent distribution implementation with direct rate normalization.
            rates = np.exp(artifact.log_rates[:, :2])
            rates /= rates.sum(axis=1, keepdims=True)
            mu = Y[:, :2].sum(axis=1)[:, None, None] * rates[None, :, :]
            reference = stats.nbinom.logpmf(Y[:, None, :2], 10., 10. / (10. + mu)).sum(axis=2)
            reference += np.log(artifact.candidate_prior)
            # Float discrepancies are reported, not judged by a convenient epsilon.
            # Independent exact rational probabilities determine the label oracle.
            expected = exact_candidates(Y, weights, artifact.candidate_prior)
            finite = np.isfinite(scores).any(axis=1)
            ties = finite & ((scores == scores.max(axis=1, keepdims=True)).sum(axis=1) > 1)
            np.testing.assert_array_equal(ties, expected == -1)
            np.testing.assert_array_equal(scores[~ties].argmax(axis=1), expected[~ties])
            changed = Y.copy()
            changed[:, 2:] = [100000, 1]
            perturbed = score_confirmation(changed, ids, artifact)
            np.testing.assert_array_equal(scores, perturbed['proposal_log_scores'])
            np.testing.assert_array_equal(result['candidate'], perturbed['candidate'])
            if args.reference:
                before = old['fixtures'][name]
                np.testing.assert_array_equal(scores, before['scores'])
                np.testing.assert_array_equal(result['candidate'][~ties], np.array(before['candidate'])[~ties])
                assert np.all(result['candidate'][ties] == -1)
                assert np.all(result['unresolved_reason'][ties] == 'tied_candidate_scores')
                assert np.all(result['unresolved_reason'][~ties] == '')
                np.testing.assert_array_equal(result['candidate'], expected)
                # Unique assignments and their contrast diagnostics must not change.
                unique = score_confirmation(Y[~ties], np.array(ids)[~ties], artifact)
                assert unique['contrasts'] == before['unique_contrasts']
            unique = score_confirmation(Y[~ties], np.array(ids)[~ties], artifact)
            records[name] = dict(scores=scores.tolist(), candidate=result['candidate'].tolist(),
                unresolved_reason=result['unresolved_reason'].tolist(), ties=int(ties.sum()),
                unique_contrasts=unique['contrasts'],
                maximum_independent_score_discrepancy=float(np.max(np.abs(scores-reference))))
        assert records['identical_proposal']['ties'] == 100
        assert records['symmetric']['ties'] == 10
        assert records['unequal_prior']['ties'] == 0
        preserved = json.loads(inputs['baseline'].read_text())
        for name, digest in preserved.items():
            if file_record(args.task/'baseline'/name)['sha256'] != digest:
                raise AssertionError('Frozen baseline changed')
        # Clone inference is outside this correction, including its retained failures.
        import haplodose
        package = Path(haplodose.__file__).parent
        for relative, digest in preserved.items():
            if relative.startswith('python/haplodose/') and relative.endswith(('.py', '.tsv')):
                if relative in ('python/haplodose/confirmation.py', 'python/haplodose/pipeline.py'):
                    continue
                assert hashlib.sha256((package/relative.removeprefix('python/haplodose/')).read_bytes()).hexdigest() == digest
        (args.out/'report.json').write_text(json.dumps(dict(fixtures=records,
            mode='candidate_audit' if args.reference else 'old_counterexample',
            old_tied_assignments=(sum(v['ties'] for v in records.values()) if not args.reference else None),
            baseline_files_preserved=len(preserved), scientific_criteria_changed=False), indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
