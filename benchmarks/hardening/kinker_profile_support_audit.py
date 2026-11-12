"""Independent scalar checks and locked-case before/after scorer regression."""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from haplodose.provenance import start_run, finish_run, file_record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    args = p.parse_args(); args.out.mkdir(exist_ok=False)
    source = Path(__file__).parents[2]
    path = source/'benchmarks/kinker/evaluate.py'
    spec = importlib.util.spec_from_file_location('support_checked_kinker', path)
    scorer = importlib.util.module_from_spec(spec); spec.loader.exec_module(scorer)
    chain = args.task/'chain_v1'
    inputs = dict(scorer=path, tests=source/'tests/test_hardening_inputs.py',
                  chain=chain/'report.json', identities=chain/'identities.json',
                  chain_manifest=chain/'run_manifest.json',
                  historical_scorer=args.task/'source_v1/benchmarks/kinker/evaluate.py')
    manifest = start_run(args.out, inputs, {},
        dict(scope='Computational scorer correction; no inference or scientific promotion'),
        mode='kinker_profile_support_audit')
    try:
        previous = json.loads(inputs['chain_manifest'].read_text())
        assert previous['status'] == 'complete'
        for group in ('input', 'output'):
            for record in previous[group]['files'].values():
                assert file_record(record['path'])['sha256'] == record['sha256']
        with (args.out/'pytest.txt').open('w') as output:
            subprocess.run([sys.executable, '-I', '-m', 'pytest', '-q', str(inputs['tests'])],
                           cwd=args.out, stdout=output, stderr=subprocess.STDOUT, check=True)
        root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
        rows = []
        for case in json.loads(inputs['chain'].read_text())['cases']:
            run = root/f"results/kinker_field_20261102_0529/haplodose_r{case['rotation']}"
            cells = pd.read_csv(run/'cells.tsv.gz', sep='\t', low_memory=False)
            query = cells.loc[cells.specimen.eq(case['case']) & cells.role.eq('query')]
            profiles = pd.read_csv(run/'profiles.tsv.gz', sep='\t')
            profiles = profiles.loc[profiles.specimen.eq(case['case'])]
            bins = pd.read_csv(run/'genes.tsv.gz', sep='\t').bin.unique()
            n = query.clone.value_counts()
            for col in ('log_ratio', 'log_dose', 'log_ratio_gene'):
                pivot = profiles.pivot(index='bin', columns='clone', values=col)
                nn = n[n.index.isin(pivot.columns)]
                old = ((pivot[nn.index]*nn.to_numpy()).sum(axis=1)/nn.sum()).reindex(bins)
                new, fraction, reason = scorer.specimen_profile(profiles, query.clone, bins, col)
                assert fraction == 1 and reason == 'complete'
                reference = np.array([math.fsum(float(pivot.loc[b, c])*int(k) for c, k in n.items())/len(query)
                                      for b in bins])
                # Standard gamma_n bound for positive-weight dot products; cancellation
                # uses the absolute sum, never relative error against a near-zero answer.
                u = np.finfo(float).eps
                operations = 2*len(n)+2
                gamma = operations*u/(1-operations*u)
                absolute = (pivot.loc[bins, n.index].abs()*n.to_numpy()).sum(axis=1).to_numpy()/len(query)
                bound = gamma*absolute
                discrepancy = np.abs(new.to_numpy()-reference)
                assert (discrepancy <= bound).all()
                assert (np.abs(old.to_numpy()-new.to_numpy()) <= 2*bound).all()
                rows.append(dict(case=case['case'], quantity=col, bins=len(bins),
                    max_scalar_discrepancy=float(discrepancy.max()),
                    max_before_after_discrepancy=float(np.abs(old.to_numpy()-new.to_numpy()).max())))
        pd.DataFrame(rows).to_csv(args.out/'regression.tsv', sep='\t', index=False)
        (args.out/'report.json').write_text(json.dumps(dict(
            complete_population_preserved=True, independent_scalar_checks=sum(r['bins'] for r in rows),
            historical_results_modified=False, no_model_or_criteria_changes=True,
            regression=rows), indent=2)+'\n')
    except BaseException as error:
        finish_run(args.out, manifest, error=error); raise
    finish_run(args.out, manifest, backend='CPU NumPy/SciPy and pytest')


if __name__ == '__main__':
    main()
