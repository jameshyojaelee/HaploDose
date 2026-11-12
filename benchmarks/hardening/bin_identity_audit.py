"""Independent named-count and no-tie replay of the bin-order correction."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.bins import make_bins
from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    prior = args.task/'probe_v1'
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        before=prior/'run_manifest.json', old_source=args.task/'probe_source_v2/source_manifest.json',
        corrected_source=args.task/'candidate_v1/source_manifest.json')
    run = start_run(args.out, inputs, {}, dict(scope='Exact input identity and computational replay only'), mode='bin_identity_audit')
    try:
        manifest = json.loads((prior/'run_manifest.json').read_text()); assert manifest['status'] == 'complete'
        for record in manifest['output']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        old_source = json.loads(inputs['old_source'].read_text())
        candidate = json.loads(inputs['corrected_source'].read_text())
        changed = sorted(name for name in old_source if candidate.get(name) != old_source[name])
        added = sorted(set(candidate)-set(old_source))
        assert changed == ['python/haplodose/bins.py']
        assert added == ['tests/test_bin_identity.py']
        for name, digest in old_source.items():
            assert file_record(args.task/'probe_source_v2'/name)['sha256'] == digest
        for name, digest in candidate.items():
            assert file_record(args.task/'candidate_v1'/name)['sha256'] == digest
        records = json.loads((prior/'cases.json').read_text())
        names = list('ABCDEF')
        counts = {'A':(2,13),'B':(3,17),'C':(5,19),'D':(101,23),'E':(7,29),'F':(11,31)}
        expected_membership = dict(A=0,B=0,C=0,D=1,E=1,F=1)
        expected_counts = np.array([[10,119],[49,83]], dtype=np.int64)
        seen = set(); new_records = []; old_variants = Counter(); old_unchanged = 0
        for record in records:
            label, order = record['condition'], list(record['input_order'])
            key = (label, ''.join(order)); assert key not in seen; seen.add(key)
            assert sorted(order) == names and label in ('unique','tied')
            starts = dict(zip(names, [1,2,3,3,4,5] if label == 'tied' else [1,2,3,4,5,6]))
            # Build annotation in an order unrelated to the input columns.
            positions = pd.DataFrame([dict(gene=g,chromosome='chr1',start=starts[g],end=starts[g]+1)
                                      for g in ['D','A','F','C','B','E']])
            bins, mapping, used = make_bins(order, positions, np.ones(6), bin_genes=3, exclude_pattern=None)
            actual = dict(zip(order, map(int, mapping)))
            assert actual == expected_membership
            assert used.gene.tolist() == names
            assert bins[['start','end','n_genes','first_gene','last_gene','genomic_block']].to_dict('list') == dict(
                start=[1,3 if label=='tied' else 4], end=[4,6 if label=='tied' else 7],
                n_genes=[3,3], first_gene=['A','D'], last_gene=['C','F'], genomic_block=[0,0])
            # Reconstruct both historical and corrected totals from ID lookup,
            # independently of the producer's column-permuted raw matrix.
            old_mapping = dict(zip(names, record['named_bin_assignment']))
            old_totals = np.array([[sum(counts[g][cell] for g in names if old_mapping[g]==b)
                                    for b in (0,1)] for cell in (0,1)])
            assert np.array_equal(old_totals, record['binned_counts'])
            corrected = np.array([[sum(counts[g][cell] for g in names if actual[g]==b)
                                   for b in (0,1)] for cell in (0,1)])
            assert np.array_equal(corrected, expected_counts)
            assert np.array_equal(corrected.sum(axis=1), old_totals.sum(axis=1))
            old_variants[label,tuple(record['named_bin_assignment'])] += 1
            unchanged = np.array_equal(corrected, old_totals)
            if label == 'unique':
                assert actual == old_mapping and unchanged
                old_unchanged += 1
            new_records.append(dict(condition=label, input_order=''.join(order), named_bin_assignment=actual,
                binned_counts=corrected.tolist(), historical_counts_identical=bool(unchanged)))
        assert len(seen) == 1440 and old_unchanged == 720
        assert sorted(n for (label,mapping),n in old_variants.items() if label=='tied') == [360,360]
        (args.out/'cases.json').write_text(json.dumps(new_records, indent=2)+'\n')
        (args.out/'report.json').write_text(json.dumps(dict(all_cases_exact=1440, count_conservation_cases=1440,
            corrected_tie_cases=720, old_unique_start_controls_exact=720,
            historical_tied_assignments=2, corrected_tied_assignments=1,
            changed_existing_source=changed, added_tests=added, old_source_files_preserved=len(old_source),
            original_tests_unchanged=True, thresholds_changed=False, scientific_accuracy_or_calibration_claim=False),indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__': main()
