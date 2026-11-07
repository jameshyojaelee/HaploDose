#!/usr/bin/env python3
"""Claim 4 SPECTRUM scRNA samples to process: eligible patients (DNA rule, results/spectrum_scorer_*/eligible_patients.tsv),
their CD45-negative and UNSORTED samples in SRP329496 (benchmarks/heldout2/spectrum_runs.tsv), kept only if the
authors' GEO cell table (GSE180661_GEO_cells.tsv.gz) labels the sample's cells. Reads metadata only.
One row per sample: sample_title, patient, group, runs (comma-separated SRR, CellRanger lanes in this order), gb,
n_geo_cells. Also writes OUT.excluded.tsv (eligible patients or samples left out, with the reason).
Usage: select_samples.py ELIGIBLE RUNS GEO_CELLS OUT
"""
import sys

import pandas as pd


def main(eligible, runs, geo, out):
    pats = pd.read_csv(eligible, sep='\t').patient
    r = pd.read_csv(runs, sep='\t')
    r['patient'] = r.sample_title.str.extract(r'^(SPECTRUM-OV-\d+)_')[0]
    r = r[r.patient.isin(pats) & r.group.isin(['CD45N', 'UNSORTED'])]
    g = pd.read_csv(geo, sep='\t', usecols=['sample']).value_counts('sample')
    s = (r.groupby(['sample_title', 'patient', 'group'])
          .agg(runs=('run', lambda x: ','.join(sorted(x))), gb=('fastq_bytes', lambda b: sum(int(v) for x in b for v in str(x).split(';')) / 1e9))
          .reset_index())
    s['n_geo_cells'] = s.sample_title.map(g).fillna(0).astype(int)
    excl = [dict(patient=p, sample_title='', reason='no CD45N or UNSORTED scRNA run in SRP329496')
            for p in pats if p not in set(s.patient)]
    excl += [dict(patient=x.patient, sample_title=x.sample_title, reason='no cells in the GEO cell table (no authors labels)')
             for x in s[s.n_geo_cells.eq(0)].itertuples()]
    s = s[s.n_geo_cells.gt(0)].sort_values(['patient', 'sample_title'])
    s.to_csv(out, sep='\t', index=False)
    pd.DataFrame(excl, columns=['patient', 'sample_title', 'reason']).to_csv(out.replace('.tsv', '') + '.excluded.tsv', sep='\t', index=False)
    print(f'{len(s)} samples, {s.patient.nunique()} patients, {s.runs.str.count(",").add(1).sum()} runs, {s.gb.sum():.0f} GB; '
          f'{len(excl)} exclusions')


if __name__ == '__main__':
    main(*sys.argv[1:5])
