"""One-fixture forensic ablation of the unchanged thin-mode recovery assertion.

Run with the already installed candidate-v1 wheel before replacing that environment. Historical
feature selection/bin boundaries are diagnostic monkeypatches confined to this process, never defaults.
No parameter tuning or frozen-result edits occur. Each ablation runs the original default simulation once.
"""
import argparse
import hashlib
import importlib.util
import inspect
import json
from pathlib import Path
import shutil
import sys
import traceback

import numpy as np
import pandas as pd

import haplodose
from haplodose import pipeline
from haplodose.clones import CloneParams
from haplodose.io import CountSource
from haplodose.simulate import write_dataset


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--baseline-source', type=Path, required=True)
    parser.add_argument('--tests-source', type=Path, required=True)
    parser.add_argument('--ablation', choices=('all', 'current', 'historical_features', 'historical_bins', 'both'),
                        default='all')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    test_path = args.tests_source / 'tests/test_split_mode.py'
    preserved = args.out / 'preserved_test_split_mode.py'
    shutil.copy2(test_path, preserved)
    assertions = load_module('preserved_thin_assertions', preserved)
    legacy_bins = load_module('haplodose._forensic_historical_bins', args.baseline_source / 'python/haplodose/bins.py')
    original_bins = pipeline.make_bins
    data = write_dataset(args.out / 'data')
    with CountSource(data / 'counts.h5') as source:
        meta = pd.read_csv(data / 'cells.tsv', sep='\t', dtype=str, keep_default_na=False)
        if 'dataset' not in meta:
            meta['dataset'] = 'dataset'
        fit = meta.role.eq('reference') & meta.ref_split.eq('fit')
        hold = meta.role.eq('reference') & meta.ref_split.eq('holdout')
        query = meta.role.isin(pipeline.QUERY_ROLES)
        types = sorted(set((meta.dataset + '|' + meta.ref_group)[fit]))
        datasets = sorted(meta.dataset.unique())
        keys = (meta.dataset + '|' + meta.ref_group).to_numpy()
        group = np.full(len(meta), -1, np.int32)
        group[fit] = pd.Index(types).get_indexer(keys[fit])
        group[query] = len(types) + pd.Index(datasets).get_indexer(meta.dataset[query])
        group[hold] = len(types) + len(datasets)
        sums = np.zeros((len(types) + len(datasets) + 1, source.n_genes))
        for a, b, X in source.chunks(97):
            sums += pipeline.grouped_column_sums(X, group[a:b], len(sums))
        historical_mean = sums.sum(axis=0) / int((query | meta.role.eq('reference')).sum())
        historical_mismatch = pipeline.type_mismatch_genes(sums, types, datasets, source.genes)
    contract = dict(counts_kind='raw', genome_build='hg38', coordinate_system='one-based-closed', gene_namespace='symbol')
    # Allows a separate isolated historical interpreter to run the same preserved assertion.
    if 'counts_kind' not in inspect.signature(pipeline.run).parameters:
        contract = {}
    variants = ('current', 'historical_features', 'historical_bins', 'both') if args.ablation == 'all' else (args.ablation,)
    records = []
    for variant in variants:
        def bins(genes, positions, mean, exclude, *positional, **kwargs):
            features = variant in ('historical_features', 'both')
            mean = historical_mean if features else mean
            exclude = set(exclude) | historical_mismatch if features else exclude
            if variant in ('historical_bins', 'both'):
                return legacy_bins.make_bins(genes, positions, mean, exclude, *positional,
                                             **{k: v for k, v in kwargs.items() if k not in
                                                ('genome_build', 'centromeres', 'max_gap_bp')})
            return original_bins(genes, positions, mean, exclude, *positional, **kwargs)
        pipeline.make_bins = bins
        destination = args.out / variant
        record = dict(ablation=variant, scientific_gate='original_thin_recovery_assertion', status='not_run')
        try:
            pipeline.run(data / 'counts.h5', data / 'cells.tsv', data / 'positions.tsv', destination,
                         params=CloneParams(split_mode='thin'), chunk=97, **contract)
            try:
                assertions.test_thin_records_mode_and_recovers_clones({'thin': destination})
                record['status'] = 'passed'
            except AssertionError:
                record.update(status='failed', traceback=traceback.format_exc())
            cells = pd.read_csv(destination / 'cells.tsv.gz', sep='\t')
            s1 = cells[cells.specimen.eq('S1') & cells.truth_clone.isin(['clone1', 'clone2']) & cells.malignant_call]
            record.update(s1_tumour_groups=int(s1.clone.nunique()), n_scored_cells=len(s1),
                          bins_sha256=digest(destination / 'bins.tsv'),
                          clone_crosstab=pd.crosstab(s1.truth_clone, s1.clone).to_dict())
        except Exception:
            record.update(status='execution_error', traceback=traceback.format_exc())
        finally:
            pipeline.make_bins = original_bins
        records.append(record)
        (args.out / 'outcomes.json').write_text(json.dumps(records, indent=2) + '\n')
    package = Path(haplodose.__file__).resolve().parent
    baseline_package = args.baseline_source / 'python/haplodose'
    hashes = {name: dict(installed=digest(package / name), baseline=digest(baseline_package / name))
              for name in ('profile.py', 'kernels/likelihood.py', 'clones.py', 'bins.py', 'pipeline.py', 'simulate.py')}
    provenance = dict(command=sys.argv, executable=sys.executable, installed_package=str(package),
                      original_assertion_sha256=digest(preserved), code_hashes=hashes,
                      likelihood_source_unchanged=hashes['kernels/likelihood.py']['installed'] ==
                          hashes['kernels/likelihood.py']['baseline'],
                      scope='one unchanged fixture; diagnostic historical toggles; no scientific claim or parameter tuning')
    (args.out / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')


if __name__ == '__main__':
    main()
