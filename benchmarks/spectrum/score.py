#!/usr/bin/env python3
"""SPECTRUM claim 4 primary endpoint (docs/STATUS.md, claim 4): per-patient clone-difference AUROC for one tool run.

Same logic as benchmarks/andor/score.py item 2, whose functions are imported, not copied:
- Truth per arm (DNA): 1 if any two DNA tumour clones differ on >= 30% of the arm's bins by >= 1 copy
  (benchmarks/andor/scdna_cn.py clone_differences).
- Tool score per arm (RNA): the largest |difference in the arm value (default median_log_ratio)| between any two of
  the tool's RNA tumour clones (andor score.py max_pair_difference). Tumour clones = clones of the query cells the
  tool calls malignant (cells.tsv.gz malignant_call), with >= max(20, 2% of the specimen's malignant query cells)
  cells. --all-cells drops the malignant filter and uses 2% of all query cells, which is exactly Andor item 2
  (cell lines have no normal cells); it exists to validate this scorer against andor/score.py.
- AUROC of the tool score for the truth (andor score.py diff_auroc).

Per specimen as in Andor. Patient (the primary unit): a patient's RNA specimens (e.g. several sites) are not pooled;
clones are compared only within a specimen, and the patient arm score is the maximum over its specimens of the
within-specimen score. Patient AUROC from that score. If no specimen has >= 2 tumour clones the score is all zero,
so the AUROC is 0.5 (the tool reports no clone difference) and the row is flagged; a missing or unreadable run is
status no_output (paired_test.py scores it 0.5, flagged). A patient whose DNA truth is all 0 or all 1 has no AUROC
(not eligible under the patient rule).

DNA input: TSV, one row per bin, columns chromosome, start, end (optional arm and bin columns are ignored and
recomputed), then one column per DNA tumour clone with its total copy number (e.g. median over the clone's cells).
BUILD (hg19 or hg38) is the build of these coordinates; arms are assigned from that build's centromeres
(p-arm end = end of the p acen band in UCSC cytoBand; hg38 uses python/haplodose/data/centromeres_hg38.tsv, which
was derived that way, hg19 uses data/external/ucsc_hg19/cytoBand.txt.gz). Chromosomes chr1-22 only; '1' and 'chr1'
both accepted. Bins with a missing value in any clone are dropped before the truth is computed.
RNA arms are the tool's own (hg38 gene positions); arm names are matched across builds.

Secondary (descriptive): RNA vs DNA tumour clone counts; precision/recall of the arms hit by the tool's own
differences.tsv calls between kept tumour clones, with the truth base rate as chance precision.

Usage: score.py TOOL_OUT DNA_CLONES PATIENT BUILD DEST [--arm-column COL] [--all-cells]
Writes DEST/scores.tsv (one row per specimen plus a level=patient row), DEST/arm_scores.tsv, DEST/run.json.
Amendment 2026-11-03 (before any SPECTRUM tool output; review F5): arm_scores.tsv gains 'measured' (>= 2 kept
clones with a value for that arm; patient row = measured in any specimen) for a descriptive common-arm
sensitivity in paired_test.py. The primary score (missing = no difference) is unchanged.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

H = Path(__file__).resolve().parents[2]
CYTOBAND = {'hg19': H / 'data/external/ucsc_hg19/cytoBand.txt.gz', 'hg38': H / 'data/external/ucsc_hg38/cytoBand.txt.gz'}
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]
COORDS = ['chromosome', 'start', 'end', 'arm', 'bin']


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


andor = _load('andor_score', H / 'benchmarks/andor/score.py')
scdna = _load('andor_scdna_cn', H / 'benchmarks/andor/scdna_cn.py')


def centromeres_from_cytoband(path):
    """p_arm_end per chromosome = end of the p-arm acen band (the rule behind centromeres_hg38.tsv)."""
    cb = pd.read_csv(path, sep='\t', header=None, names=['chromosome', 'start', 'end', 'band', 'stain'])
    p = cb[cb.stain.eq('acen') & cb.band.str.startswith('p')]
    return p.groupby('chromosome').end.max().rename('p_arm_end').reset_index()


def centromere_file(build, dest):
    if build == 'hg38':
        return andor.CEN
    path = Path(dest) / f'centromeres_{build}.tsv'
    centromeres_from_cytoband(CYTOBAND[build]).to_csv(path, sep='\t', index=False)
    return path


def dna_truth(dna_path, build, dest):
    """DNA tumour clone names, arm order, truth arm set and the difference rows."""
    d = pd.read_csv(dna_path, sep='\t')
    names = [c for c in d.columns if c not in COORDS]
    d['chromosome'] = d.chromosome.astype(str).map(lambda c: c if c.startswith('chr') else f'chr{c}')
    d = d[d.chromosome.isin(AUTOSOMES)]
    n_bins_in = len(d)
    d = d.dropna(subset=names).reset_index(drop=True)
    arm_id = andor.arm_of(d.chromosome, d.start, centromere_file(build, dest))
    diffs = scdna.clone_differences(d[names], names, arm_id)
    return names, pd.Index(pd.unique(arm_id)), set(diffs.arm), diffs, dict(dna_bins=n_bins_in, dna_bins_used=len(d))


def is_true(x):
    return x.astype(str).str.lower().eq('true')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('tool_out')
    ap.add_argument('dna_clones')
    ap.add_argument('patient')
    ap.add_argument('build', choices=['hg19', 'hg38'])
    ap.add_argument('dest')
    ap.add_argument('--arm-column', default='median_log_ratio')
    ap.add_argument('--all-cells', action='store_true', help='no malignant filter (cell lines; Andor validation)')
    a = ap.parse_args(argv)
    tool_out, dest = Path(a.tool_out), Path(a.dest)
    if dest.exists() and any(dest.iterdir()):
        raise FileExistsError(dest)
    dest.mkdir(parents=True, exist_ok=True)
    names, all_arms, diff_arms, dna_diffs, dna_qc = dna_truth(a.dna_clones, a.build, dest)
    dna_diffs.to_csv(dest / 'dna_clone_differences.tsv', sep='\t', index=False)
    truth = pd.Series(all_arms.isin(diff_arms).astype(int), index=all_arms)
    base = dict(patient=a.patient, arm_column=a.arm_column, all_cells=a.all_cells, dna_clones=len(names),
                dna_subclonal_arms=int(truth.sum()), arms=len(all_arms), **dna_qc)
    run = dict(tool_out=str(tool_out.resolve()), dna_table=str(Path(a.dna_clones).resolve()), build=a.build, **base)

    need = ['cells.tsv.gz', 'arms.tsv']
    if not all((tool_out / f).exists() for f in need):
        pd.DataFrame([dict(level='patient', specimen=None, status='no_output', diff_auroc=None, **base)]).to_csv(
            dest / 'scores.tsv', sep='\t', index=False)
        (dest / 'run.json').write_text(json.dumps(dict(run, status='no_output', missing=need), indent=2) + '\n')
        print('no_output', tool_out)
        return

    cells = pd.read_csv(tool_out / 'cells.tsv.gz', sep='\t', dtype={'clone': str})
    arms = pd.read_csv(tool_out / 'arms.tsv', sep='\t')
    dpath = tool_out / 'differences.tsv'
    calls = pd.read_csv(dpath, sep='\t') if dpath.exists() and dpath.read_text().strip() else pd.DataFrame()
    q = cells[cells.role.eq('query') & cells.patient.astype(str).eq(a.patient)]
    if not a.all_cells:
        q = q[is_true(q.malignant_call)]
    rows, arm_rows = [], []
    patient_score = pd.Series(0.0, index=all_arms)
    patient_measured = pd.Series(False, index=all_arms)
    called_all, n_max = set(), 0
    for spec, c in q.groupby('specimen'):
        sizes = c.clone[c.clone.notna() & (c.clone != 'NA')].value_counts()
        keep = sizes[sizes >= max(20, 0.02 * len(c))]
        names_rna = list(keep.index)
        piv = arms[arms.specimen.eq(spec) & arms.clone.isin(names_rna)].pivot(index='arm', columns='clone',
                                                                             values=a.arm_column).reindex(all_arms)
        score = andor.max_pair_difference(piv, names_rna, all_arms)
        measured = (piv[names_rna].notna().sum(axis=1) >= 2) if len(names_rna) else pd.Series(False, index=all_arms)
        patient_measured |= measured
        patient_score = np.maximum(patient_score, score)
        n_max = max(n_max, len(names_rna))
        called = set()
        if len(calls):
            dd = calls[calls.specimen.eq(spec) & calls.clone_a.isin(names_rna) & calls.clone_b.isin(names_rna)]
            for r in dd.itertuples():
                called |= set(andor.arm_of([r.chromosome, r.chromosome], [r.start, r.end]))
        called &= set(all_arms)
        called_all |= called
        tp = len(called & diff_arms)
        rows.append(dict(level='specimen', specimen=spec, status='ok', n_cells=len(c), rna_tumour_clones=len(names_rna),
                         diff_auroc=andor.diff_auroc(truth, score, len(names_rna)), diff_called_arms=len(called),
                         diff_called_true=tp, **base))
        arm_rows.append(pd.DataFrame({'specimen': spec, 'arm': all_arms, 'truth': truth.to_numpy(), 'score': score.to_numpy(),
                                     'measured': measured.to_numpy()}))
    auroc = andor.diff_auroc(truth, patient_score, n_max)
    flag = None
    if auroc is None and 0 < truth.sum() < len(truth):     # < 2 tumour clones in every specimen: no difference reported
        auroc, flag = 0.5, 'fewer_than_2_tumour_clones'
    tp = len(called_all & diff_arms)
    rows.append(dict(level='patient', specimen=None, status='ok', n_cells=len(q), rna_tumour_clones=n_max,
                     diff_auroc=auroc, flag=flag, diff_called_arms=len(called_all), diff_called_true=tp,
                     diff_precision=tp / len(called_all) if called_all else None,
                     diff_recall=tp / len(diff_arms) if diff_arms else None,
                     diff_chance_precision=len(diff_arms) / len(all_arms), **base))
    arm_rows.append(pd.DataFrame({'specimen': 'patient', 'arm': all_arms, 'truth': truth.to_numpy(), 'score': patient_score.to_numpy(),
                                 'measured': patient_measured.to_numpy()}))
    t = pd.DataFrame(rows)
    t.to_csv(dest / 'scores.tsv', sep='\t', index=False)
    pd.concat(arm_rows).to_csv(dest / 'arm_scores.tsv', sep='\t', index=False)
    (dest / 'run.json').write_text(json.dumps(dict(run, status='ok'), indent=2) + '\n')
    print(t.to_string())


if __name__ == '__main__':
    main()
