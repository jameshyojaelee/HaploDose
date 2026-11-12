"""Exhaustive feature-order challenge of tied-start bin membership and raw totals."""
import argparse
from collections import Counter
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.bins import make_bins
from haplodose.provenance import start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    run = start_run(args.out, dict(script=Path(__file__), design=args.task/'DESIGN.md'),
        dict(permutations=720), dict(scope='Exact computational identity, not scientific accuracy'), mode='bin_identity_probe')
    try:
        names = np.array(list('ABCDEF'))
        counts = np.array([[2,3,5,101,7,11],[13,17,19,23,29,31]], dtype=np.int64)
        results = {}; rows = []
        for label, starts in [('tied', [1,2,3,3,4,5]), ('unique', [1,2,3,4,5,6])]:
            positions = pd.DataFrame(dict(gene=names, chromosome='chr1', start=starts, end=np.array(starts)+1))
            variants = Counter(); total_checks = 0
            for order in itertools.permutations(range(6)):
                perm = np.array(order)
                bins, mapping, used = make_bins(names[perm], positions.iloc[::-1], np.ones(6),
                    bin_genes=3, exclude_pattern=None)
                by_name = dict(zip(names[perm], map(int, mapping)))
                allocated = np.zeros((2, len(bins)), dtype=np.int64)
                for j, k in enumerate(mapping):
                    allocated[:, k] += counts[:, perm[j]]
                assert np.array_equal(allocated.sum(axis=1), counts.sum(axis=1))
                assert dict(zip(used.gene, map(int, used.bin))) == by_name
                key = tuple(by_name[name] for name in names); variants[key] += 1; total_checks += 1
                rows.append(dict(condition=label, input_order=''.join(names[perm]), named_bin_assignment=list(key),
                    binned_counts=allocated.tolist()))
            expected = (0,0,0,1,1,1)
            results[label] = dict(permutations=total_checks, unique_named_assignments=len(variants),
                canonical_reference_mapping=list(expected), canonical_matches=variants[expected],
                variants=[dict(mapping=list(key), count=n) for key,n in sorted(variants.items())],
                exact_count_conservation=True)
        assert results['unique']['canonical_matches'] == 720
        (args.out/'cases.json').write_text(json.dumps(rows, indent=2)+'\n')
        (args.out/'report.json').write_text(json.dumps(dict(conditions=results, scientific_threshold_changed=False,
            biological_validation=False), indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__': main()
