"""Replay both bulk modes with the preceding wheel; require exact numeric parity."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import haplodose
from haplodose.bulk import run_bulk
from haplodose.provenance import start_run, finish_run, file_record


def verify_artifact(path):
    artifact = json.loads(path.read_text())
    import haplodose._core as native
    if file_record(artifact['wheel']['path'])['sha256'] != artifact['wheel']['sha256']:
        raise AssertionError('Changed wheel')
    if file_record(native.__file__)['sha256'] != artifact['native']['sha256']:
        raise AssertionError('Unexpected native binary')
    if Path(haplodose.__file__).resolve() != Path(artifact['package']['path']).resolve():
        raise AssertionError('Unexpected package import')
    return artifact


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--validation', type=Path, required=True)
    parser.add_argument('--old-artifact', type=Path, required=True)
    parser.add_argument('--mode', choices=('audit', 'reference'), default='audit')
    parser.add_argument('--case', type=Path)
    args = parser.parse_args()
    if args.mode == 'reference':
        verify_artifact(args.old_artifact)
        baseline = json.loads((args.task/'baseline/source_manifest.json').read_text())
        package = Path(haplodose.__file__).parent
        for name, digest in baseline.items():
            if name.startswith('python/haplodose/') and name.endswith(('.py', '.tsv')):
                if file_record(package/name.removeprefix('python/haplodose/'))['sha256'] != digest:
                    raise AssertionError('Reference package differs from pre-edit baseline')
        configuration = json.loads((args.case/'out/run.json').read_text())
        run_bulk(args.case/'counts.tsv', args.case/'samples.tsv', args.case/'positions.tsv',
            args.out, bin_stat=configuration['bin_stat'], counts_kind='raw', genome_build='hg38',
            coordinate_system='zero-based-half-open', gene_namespace='custom')
        return
    verify_artifact(args.validation/'installed_artifact.json')
    args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        candidate_artifact=args.validation/'installed_artifact.json', old_artifact=args.old_artifact,
        candidate_pytest=args.validation/'pytest.txt', baseline=args.task/'baseline/source_manifest.json')
    # pytest also leaves a "current" symlink to the last case directory.
    cases = sorted({p.resolve() for p in (args.validation/'pytest_tmp').glob(
        'test_bulk_recovers_arm_events_*/out/run.json')})
    if len(cases) != 2:
        raise AssertionError('Expected both executed bulk-mode fixtures')
    for index, configuration in enumerate(cases):
        inputs[f'{index}:configuration'] = configuration
        inputs[f'{index}:manifest'] = configuration.parent/'run_manifest.json'
    record = start_run(args.out, inputs, dict(comparison='exact; no numeric tolerance'),
        dict(scope='Output-unit correction only; original scientific assertions unchanged'),
        mode='bulk_units_parity')
    try:
        reports = []
        modes = set()
        old = json.loads(args.old_artifact.read_text())
        for configuration in cases:
            case = configuration.parent.parent
            mode = json.loads(configuration.read_text())['bin_stat']
            if mode in modes:
                raise AssertionError('Duplicate bulk mode')
            modes.add(mode)
            target = args.out/('reference_' + mode)
            command = [old['python'], '-I', str(Path(__file__).resolve()), '--mode', 'reference',
                '--task', str(args.task), '--validation', str(args.validation), '--old-artifact',
                str(args.old_artifact), '--case', str(case), '--out', str(target)]
            with (args.out/(mode + '_reference.log')).open('w') as stream:
                subprocess.run(command, check=True, stdout=stream, stderr=subprocess.STDOUT)
            pairs = []
            for name in ('arms.tsv', 'segments.tsv', 'samples.tsv', 'bins.tsv', 'genes.tsv.gz'):
                before = pd.read_csv(target/name, sep='\t')
                after = pd.read_csv(case/'out'/name, sep='\t')
                if name == 'segments.tsv':
                    units = 'median_derived_pseudocounts' if mode == 'median' else 'input_count_units'
                    if not before.count_units.eq('input_count_units').all() or not after.count_units.eq(units).all():
                        raise AssertionError('Incorrect count units')
                    before = before.drop(columns='count_units')
                    after = after.drop(columns='count_units')
                if name == 'samples.tsv':
                    added = {'binned_count_units', 'input_binned_counts', 'input_binned_count_units'}
                    if set(after) - set(before) != added:
                        raise AssertionError('Unexpected sample-schema difference')
                    genes = pd.read_csv(target/'genes.tsv.gz', sep='\t').gene
                    counts = pd.read_csv(case/'counts.tsv', sep='\t', index_col=0)
                    # Python integer summation, separate from production bin aggregation.
                    total = [sum(int(v) for v in counts.loc[genes, sample]) for sample in after['sample']]
                    np.testing.assert_array_equal(after.input_binned_counts, total)
                    after = after.drop(columns=list(added))
                pd.testing.assert_frame_equal(before, after, check_exact=True)
                pairs.append(dict(file=name, before=file_record(target/name), after=file_record(case/'out'/name)))
            np.testing.assert_array_equal(np.load(target/'log_ratio.npy'), np.load(case/'out/log_ratio.npy'))
            before = json.loads((target/'run.json').read_text())
            after = json.loads(configuration.read_text())
            for key in ('binned_count_units', 'input_binned_count_units'):
                after.pop(key)
            if before != after:
                raise AssertionError('Unexpected run metadata difference')
            reports.append(dict(mode=mode, exact_original_tables=pairs, exact_profile_array=True,
                reference_manifest=file_record(target/'run_manifest.json')))
        if modes != {'sum', 'median'}:
            raise AssertionError('Missing bulk mode')
        baseline = json.loads(inputs['baseline'].read_text())
        for name, digest in baseline.items():
            if hashlib.sha256((args.task/'baseline'/name).read_bytes()).hexdigest() != digest:
                raise AssertionError('Frozen baseline changed')
        (args.out/'report.json').write_text(json.dumps(dict(cases=reports,
            baseline_files=len(baseline), scientific_criteria_changed=False), indent=2))
        finish_run(args.out, record, backend='cpu')
    except Exception as error:
        finish_run(args.out, record, error=error)
        raise


if __name__ == '__main__':
    main()
