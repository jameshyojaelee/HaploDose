"""Corrected development-DNA truth and depth diagnostics; run on a compute node.

Writes only a fresh --out directory. Reads existing Andor binned DNA and small
SPECTRUM metadata; no downloads or historical-result changes. Depth repetitions
are simulations of one line, not independent biological replication.
"""
import argparse
import importlib.util
import json
import platform
import shutil
import sys
import traceback
from pathlib import Path

import pandas as pd

H = Path(__file__).resolve().parents[2]


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, H / relative)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--data-root', type=Path, default=H,
                    help='live repository with data/results/dev sources; H remains the frozen code root')
    ap.add_argument('--mode', choices=['metadata', 'andor', 'all'], default='all')
    ap.add_argument('--depth-repeats', type=int, default=1)
    ap.add_argument('--seed', type=int, default=20261106)
    a = ap.parse_args(argv)
    data_root = a.data_root.resolve()
    if not 1 <= a.depth_repeats <= 10:
        raise ValueError('depth repeats must be between 1 and 10')
    out = a.out.resolve()
    if out.is_relative_to(data_root / 'data') or out.is_relative_to(data_root / 'results'):
        raise ValueError('use a fresh diagnostic root outside data/ and historical results/')
    out.mkdir(parents=True, exist_ok=False)
    stages = []

    def stage(name, action):
        try:
            action()
            stages.append(dict(stage=name, status='ok'))
        except Exception as exc:
            stages.append(dict(stage=name, status='failed', error=str(exc), traceback=traceback.format_exc()))
        (out / 'stages.json').write_text(json.dumps(stages, indent=2) + '\n')

    if a.mode in ['metadata', 'all']:
        corr = load('hardening_correspondence', 'benchmarks/spectrum/correspondence.py')
        stage('spectrum_metadata', lambda: corr.main([
            str(data_root / 'dev/spectrum_run_plan_20261103_0750/samples.tsv'),
            str(data_root / 'data/external/spectrum_ov_scdna_synapse'), str(out / 'spectrum_metadata')]))
    if a.mode in ['andor', 'all']:
        full = load('hardening_full_dna', 'benchmarks/andor/scdna_cn.py')
        low = load('hardening_low_dna', 'benchmarks/andor/scdna_lowdepth.py')
        sim = load('hardening_sim_dna', 'benchmarks/andor/simulate_lowdepth.py')
        compare = load('hardening_compare_dna', 'benchmarks/andor/compare_depth.py')
        source = data_root / 'data/external/andor2020_scdna'
        baseline = data_root / 'results/andor_20261101_1224'
        stage('andor_full_SNU-601', lambda: full.main([
            '--bins-dir', str(source / 'bins_original_bam'), '--lines', 'SNU-601', '--out', str(out / 'andor_full')]))
        for line in ['MKN-45', 'NUGC-4']:
            stage(f'andor_low_{line}', lambda line=line: low.main([
                '--bins-dir', str(source / 'bins'), '--lines', line, '--out', str(out / f'andor_low_{line}')]))
            stage(f'andor_old_new_{line}', lambda line=line: compare.main(
                baseline / 'scdna_lowdepth' / line, out / f'andor_low_{line}' / line, out / f'old_new_{line}'))
        stage('andor_old_new_SNU-601', lambda: compare.main(
            baseline / 'scdna_fulldepth/SNU-601', out / 'andor_full/SNU-601', out / 'old_new_SNU-601'))
        for repeat in range(a.depth_repeats):
            root = out / f'depth_{repeat}'
            root.mkdir()
            shutil.copy2(source / 'bin_annotation_500kb.tsv', root / 'bin_annotation_500kb.tsv')
            sim.SEED = a.seed + repeat
            stage(f'depth_simulate_{repeat}', lambda root=root: sim.main(
                source / 'bins_original_bam', source / 'bins', 'SNU-601', root / 'sim_bins'))
            stage(f'depth_infer_{repeat}', lambda root=root: low.main([
                '--bins-dir', str(root / 'sim_bins'), '--lines', 'SNU-601', '--out', str(root / 'inferred')]))
            stage(f'depth_compare_{repeat}', lambda root=root: compare.main(
                out / 'andor_full/SNU-601', root / 'inferred/SNU-601', root / 'comparison'))
        rows = []
        for path in sorted(out.glob('**/qc.json')):
            qc = json.loads(path.read_text())
            rows.append(dict(path=str(path.relative_to(out)), line=qc['line'],
                             dna_cells=qc['n_pass_mapd'], dna_clones=len(qc['clones']),
                             median_fragments=qc['median_fragments'], callable_bp=qc.get('callable_bp'),
                             uncertainty=qc.get('truth_uncertainty')))
        pd.DataFrame(rows).to_csv(out / 'andor_truth_trace.tsv', sep='\t', index=False)
    summary = dict(schema_version=2, args={k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()},
                   python=sys.version, host=platform.node(), stages=stages,
                   biological_units='3 development lines; depth repetitions are not biological replication',
                   passed=all(x['status'] == 'ok' for x in stages),
                   limits='Andor RNA/DNA cells are separate. Low-depth single-cluster results do not certify clonality. '
                          'SPECTRUM metadata cannot establish same-biopsy matching.')
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return 0 if summary['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
