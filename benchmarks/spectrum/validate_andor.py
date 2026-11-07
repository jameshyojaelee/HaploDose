#!/usr/bin/env python3
"""Validate the SPECTRUM claim 4 scorer on Andor 2020 SNU-601 development data before any SPECTRUM output exists.

Checks (each printed PASS/FAIL; exit code 1 if any fails):
1. Centromeres: the cytoBand rule in score.py (end of the p acen band) reproduces
   python/haplodose/data/centromeres_hg38.tsv from data/external/ucsc_hg38/cytoBand.txt.gz.
2. DNA truth: score.py dna_truth on SNU-601 clone_bins.tsv (scdna_combined) recomputes the stored arm column and
   reproduces the stored clone_differences.tsv row for row (this also checks the clone_differences refactor of
   benchmarks/andor/scdna_cn.py).
3. AUROC match: for each SNU-601 run, andor/score.py (refactored) reproduces its stored diff_auroc, and spectrum
   score.py --all-cells gives the same diff_auroc, specimen row and patient row (SNU-601 has one specimen), with
   exact float equality. The malignant-filter mode (the SPECTRUM primary mode) is run too and reported, not checked
   (a cell line has no normal cells, so its value need not match).
4. Synthetic cases with known answers: malignant filter and size threshold, < 2 tumour clones -> 0.5 flagged,
   patient score = max over specimens, hg19 arms, no_output; paired_test.py on a toy score tree.

Usage: validate_andor.py OUT
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

H = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(H / 'benchmarks/spectrum'))
import score as sp          # noqa: E402
import paired_test as pt    # noqa: E402

A = H / 'results/andor_20261101_1224'
SCDNA = A / 'scdna_combined'
LINE = 'SNU-601'
RUNS = {  # tool -> (run dir, stored andor score.py output)
    'haplodose_frozen1424': (H / 'results/andor_frozen1424_20261102_1542/v6_s20261031/runs/SNU-601/haplodose',
                             H / 'results/andor_frozen1424_20261102_1542/v6_s20261031/scores/SNU-601/haplodose/scores.tsv'),
    'haplodose_1101': (A / 'rna_runs/SNU-601/haplodose', A / 'scores_combined/SNU-601/haplodose/scores.tsv'),
    'infercnvpy': (A / 'rna_runs/SNU-601/infercnvpy_20261101_122809', A / 'scores_combined/SNU-601/infercnvpy/scores.tsv'),
    'infercnv_r': (A / 'rna_runs/SNU-601/infercnv_r_20261101_122851', A / 'scores_combined/SNU-601/infercnv_r/scores.tsv'),
    'prototype': (A / 'rna_runs/SNU-601/prototype_20261101_122802', A / 'scores_combined/SNU-601/prototype/scores.tsv'),
}
results = []


def check(name, ok, detail=''):
    results.append(dict(check=name, result='PASS' if ok else 'FAIL', detail=str(detail)))
    print(('PASS ' if ok else 'FAIL ') + name, detail)


def patient_row(dest):
    t = pd.read_csv(Path(dest) / 'scores.tsv', sep='\t')
    return t[t.level.eq('patient')].iloc[0], t[t.level.eq('specimen')]


def synthetic(out):
    """Toy tool run and DNA table with known answers."""
    rng = np.random.default_rng(1)
    arms_list = ['chr1p', 'chr1q', 'chr2p', 'chr2q', 'chr3p', 'chr3q']
    # hg19 chr1 p-arm end 125,000,000; hg38 123,400,000. A bin at 124,000,000 is chr1q in hg38, chr1p in hg19.
    cen19 = sp.centromeres_from_cytoband(sp.CYTOBAND['hg19']).set_index('chromosome').p_arm_end
    bins = []
    for ch in ['1', '2', '3']:
        e = int(cen19[f'chr{ch}'])
        bins += [(ch, s, s + 500_000) for s in range(e - 10_000_000, e - 5_000_000, 500_000)]   # p arm
        bins += [(ch, s, s + 500_000) for s in range(e + 5_000_000, e + 10_000_000, 500_000)]   # q arm
    dna = pd.DataFrame(bins, columns=['chromosome', 'start', 'end'])
    dna['d1'] = 2.0
    dna['d2'] = 2.0
    dna.loc[(dna.chromosome == '2') & (dna.start > cen19['chr2']), 'd2'] = 3.0      # chr2q differs on all bins
    dna.loc[dna.index[(dna.chromosome == '3') & (dna.start < cen19['chr3'])][:2], 'd2'] = 1.0   # chr3p 2/10 bins: below 30%
    dna['d3'] = dna.d1
    dna.loc[2, 'd3'] = np.nan                                                              # dropped bin
    dna_path = out / 'toy_dna.tsv'
    dna.to_csv(dna_path, sep='\t', index=False)
    names, all_arms, diff_arms, _, qc = sp.dna_truth(dna_path, 'hg19', out)
    check('synthetic: hg19 arms and 30% rule', list(all_arms) == arms_list and diff_arms == {'chr2q'} and qc['dna_bins_used'] == len(dna) - 1,
          dict(arms=list(all_arms), diff=sorted(diff_arms)))

    def tool_run(d, layout):
        d.mkdir(parents=True)
        cells, arms = [], []
        for spec, clones in layout.items():
            for cl, (n, malig, values) in clones.items():
                cells += [dict(cell_id=f'{spec}_{cl}_{i}', specimen=spec, patient='P1', role='query', clone=cl,
                               malignant_call=malig) for i in range(n)]
                arms += [dict(specimen=spec, clone=cl, arm=a_, median_log_ratio=v) for a_, v in zip(arms_list, values)]
        cells.append(dict(cell_id='ref1', specimen='REF', patient='P1', role='reference', clone='NA', malignant_call=False))
        pd.DataFrame(cells).to_csv(d / 'cells.tsv.gz', sep='\t', index=False, compression='gzip')
        pd.DataFrame(arms).to_csv(d / 'arms.tsv', sep='\t', index=False)
        (d / 'differences.tsv').write_text('\n')

    flat = [0.0] * 6
    gain2q = [0.0, 0.0, 0.0, 0.4, 0.0, 0.0]
    big_normal = [0.0, 0.0, 0.0, 0.0, 0.9, 0.0]   # a normal clone that differs at chr3p: must be ignored
    # Specimen S1: tumour clones a (100) and b (50) differ at chr2q; normal clone n (200, malignant False).
    # Specimen S2: tumour clones c (100), d (5 cells: below max(20, 2%)) differ at chr1p; ignored.
    tool_run(out / 'toy_run', {'S1': {'a': (100, True, flat), 'b': (50, True, gain2q), 'n': (200, False, big_normal)},
                               'S2': {'c': (100, True, flat), 'd': (5, True, [0.9] + [0.0] * 5)}})
    sp.main(['--arm-column', 'median_log_ratio', str(out / 'toy_run'), str(dna_path), 'P1', 'hg19', str(out / 'toy_score')])
    pr, specs = patient_row(out / 'toy_score')
    arm_s = pd.read_csv(out / 'toy_score/arm_scores.tsv', sep='\t').query('specimen == "patient"').set_index('arm').score
    check('synthetic: malignant filter, size threshold, patient max over specimens',
          pr.diff_auroc == 1.0 and pr.rna_tumour_clones == 2 and abs(arm_s['chr2q'] - 0.4) < 1e-12 and arm_s.drop('chr2q').eq(0).all()
          and specs.set_index('specimen').rna_tumour_clones.to_dict() == {'S1': 2, 'S2': 1},
          dict(auroc=pr.diff_auroc, score=arm_s.round(3).to_dict()))
    sp.main([str(out / 'toy_run'), str(dna_path), 'P1', 'hg19', str(out / 'toy_score_allcells'), '--all-cells'])
    pr2, _ = patient_row(out / 'toy_score_allcells')
    check('synthetic: --all-cells keeps the normal clone (score at chr3p rises, AUROC falls)', pr2.diff_auroc < 1.0,
          dict(auroc=pr2.diff_auroc))
    tool_run(out / 'toy_one_clone', {'S1': {'a': (100, True, flat), 'n': (200, False, big_normal)}})
    sp.main([str(out / 'toy_one_clone'), str(dna_path), 'P1', 'hg19', str(out / 'toy_score_one')])
    pr3, _ = patient_row(out / 'toy_score_one')
    check('synthetic: < 2 tumour clones -> AUROC 0.5 flagged', pr3.diff_auroc == 0.5 and pr3.flag == 'fewer_than_2_tumour_clones',
          dict(auroc=pr3.diff_auroc, flag=pr3.flag))
    sp.main([str(out / 'toy_missing'), str(dna_path), 'P1', 'hg19', str(out / 'toy_score_missing')])
    pr4, _ = patient_row(out / 'toy_score_missing')
    check('synthetic: missing run -> status no_output', pr4.status == 'no_output', dict(status=pr4.status))

    # paired_test: 8 patients; haplodose 0.8, numbat 0.6 - 0.01 i (distinct differences, so the Wilcoxon is exact); infercnvpy equal to haplodose; copykat no output;
    # scevan missing for 2 patients (scored 0.5).
    root = out / 'toy_tree'
    for i in range(8):
        for tool, v in [('haplodose', 0.8), ('numbat', 0.6 - 0.01 * i), ('infercnvpy', 0.8), ('scevan', 0.7)]:
            if tool == 'scevan' and i < 2:
                continue
            d = root / f'P{i}' / tool
            d.mkdir(parents=True)
            pd.DataFrame([dict(level='patient', status='ok', diff_auroc=v)]).to_csv(d / 'scores.tsv', sep='\t', index=False)
    pd.DataFrame({'patient': [f'P{i}' for i in range(8)]}).to_csv(out / 'toy_patients.tsv', sep='\t', index=False)
    pt.main(out / 'toy_patients.tsv', root, out / 'toy_paired')
    s = json.loads((out / 'toy_paired/primary.json').read_text())
    t = {x['comparator']: x for x in s['tests']}
    exp_scevan = np.mean([0.3, 0.3] + [0.1] * 6)
    check('synthetic: paired_test family, exact Wilcoxon, CI, missing = 0.5',
          s['family'] == ['infercnvpy', 'copykat', 'scevan', 'numbat'] and s['no_output_any_patient'] == ['copykat']
          and s['completion_ok'] == dict(infercnvpy=True, copykat=False, scevan=True, numbat=True)
          and t['copykat']['superior'] and s['sensitivity_both_scored'][1]['n'] == 0
          and abs(t['numbat']['mean_diff'] - 0.235) < 1e-12 and abs(t['numbat']['p'] - 2 / 256) < 1e-12
          and np.isnan(t['infercnvpy']['p']) and t['infercnvpy']['p_holm'] == 1.0 and not t['infercnvpy']['superior']
          and abs(t['scevan']['mean_diff'] - exp_scevan) < 1e-12 and s['failures_scored_0_5']['scevan'] == 2
          and t['numbat']['superior'] and not s['sota'],
          dict(family=s['family'], numbat_p=t['numbat']['p'], numbat_holm=t['numbat']['p_holm'],
               scevan=t['scevan']['mean_diff'], sota=s['sota']))

    # Review F5: score.py marks measured arms; S2 has 1 kept tumour clone, so only S1 arms are measured.
    m = pd.read_csv(out / 'toy_score/arm_scores.tsv', sep='\t').query('specimen == "patient"').set_index('arm').measured
    check('synthetic: arm_scores measured column (>= 2 kept clones with a value)', bool(m.all()), dict(measured=m.to_dict()))
    # Common-arm AUROC: truth 1,1,0,0,0; comparator missing a3,a4 (filled 0 in the primary, which would favour it).
    croot = out / 'toy_common'
    arms5 = ['a1', 'a2', 'a3', 'a4', 'a5']
    for tool, sc, ms in [('haplodose', [1, 1, 2, 2, 0], [1] * 5), ('infercnvpy', [1, 1, 0, 0, .5], [1, 1, 0, 0, 1])]:
        d = croot / 'P0' / tool
        d.mkdir(parents=True)
        pd.DataFrame(dict(specimen='patient', arm=arms5, truth=[1, 1, 0, 0, 0], score=sc, measured=[bool(x) for x in ms])
                     ).to_csv(d / 'arm_scores.tsv', sep='\t', index=False)
    ca = pt.common_arm_auroc(croot, 'P0', 'haplodose', 'infercnvpy')
    check('synthetic: common-arm AUROC uses arms both tools measured', ca == (1.0, 1.0), dict(common=ca))


def main(out):
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(out)
    out.mkdir(parents=True, exist_ok=True)

    cen = pd.read_csv(sp.andor.CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    der = sp.centromeres_from_cytoband(sp.CYTOBAND['hg38']).set_index('chromosome').p_arm_end.reindex(cen.index)
    check('1 centromeres: cytoBand rule reproduces centromeres_hg38.tsv', bool((der == cen).all()),
          f'{int((der == cen).sum())}/{len(cen)} chromosomes equal')

    cb_path = SCDNA / LINE / 'clone_bins.tsv'
    names, all_arms, diff_arms, diffs, _ = sp.dna_truth(cb_path, 'hg38', out)
    cb = pd.read_csv(cb_path, sep='\t')
    stored = pd.read_csv(SCDNA / LINE / 'clone_differences.tsv', sep='\t')
    arm_recomputed = sp.andor.arm_of(cb.chromosome, cb.start)
    check('2a arms recomputed from coordinates equal the stored arm column', bool((arm_recomputed == cb.arm.to_numpy()).all()))
    x, y = diffs.reset_index(drop=True), stored.reset_index(drop=True)
    # Labels exact; floats to 1e-12 (the stored file went through to_csv/read_csv, which can move the last ulp).
    same = (len(x) == len(y) and all((x[c] == y[c]).all() for c in ['clone_a', 'clone_b', 'arm'])
            and all(np.allclose(x[c], y[c], rtol=0, atol=1e-12) for c in ['frac_bins_diff', 'mean_diff'])
            and diff_arms == set(y.arm))
    check('2b clone_differences (refactored) reproduces stored clone_differences.tsv (rows, labels, truth arm set)', same,
          f'{len(diffs)} rows vs {len(stored)}; truth arms {len(diff_arms)}/{len(all_arms)}')
    qc = json.loads((SCDNA / LINE / 'qc.json').read_text())
    check('2c DNA clones in clone_bins.tsv = qc.json clones', sorted(names) == sorted(qc['clones']), len(names))

    rows = []
    for tool, (run, stored_scores) in RUNS.items():
        old = pd.read_csv(stored_scores, sep='\t').iloc[0]
        d_andor, d_all, d_mal = out / 'andor' / tool, out / 'spectrum_all_cells' / tool, out / 'spectrum_malignant' / tool
        sp.andor.main(str(run), str(SCDNA), LINE, str(d_andor))
        new = pd.read_csv(d_andor / 'scores.tsv', sep='\t').iloc[0]
        sp.main([str(run), str(cb_path), LINE, 'hg38', str(d_all), '--all-cells'])
        pa, sa = patient_row(d_all)
        sp.main([str(run), str(cb_path), LINE, 'hg38', str(d_mal)])
        pm, _ = patient_row(d_mal)
        rows.append(dict(tool=tool, stored_andor=old.diff_auroc, andor_refactored=new.diff_auroc,
                         spectrum_all_cells_specimen=sa.diff_auroc.iloc[0], spectrum_all_cells_patient=pa.diff_auroc,
                         spectrum_malignant_patient=pm.diff_auroc, rna_clones_andor=int(old.rna_clones),
                         rna_tumour_clones_all=int(pa.rna_tumour_clones), rna_tumour_clones_malignant=int(pm.rna_tumour_clones)))
        r = rows[-1]
        check(f'3 {tool}: andor refactor = stored', r['andor_refactored'] == r['stored_andor'], r['stored_andor'])
        check(f'3 {tool}: spectrum --all-cells specimen and patient = stored andor item 2',
              r['spectrum_all_cells_specimen'] == r['stored_andor'] == r['spectrum_all_cells_patient']
              and r['rna_tumour_clones_all'] == r['rna_clones_andor'],
              (r['spectrum_all_cells_specimen'], r['spectrum_all_cells_patient']))
    pd.DataFrame(rows).to_csv(out / 'andor_match.tsv', sep='\t', index=False)
    print(pd.DataFrame(rows).to_string())

    (out / 'synthetic').mkdir()
    synthetic(out / 'synthetic')
    pd.DataFrame(results).to_csv(out / 'checks.tsv', sep='\t', index=False)
    n_fail = sum(r['result'] == 'FAIL' for r in results)
    print(f'{len(results) - n_fail}/{len(results)} checks passed')
    sys.exit(1 if n_fail else 0)


if __name__ == '__main__':
    main(sys.argv[1])
