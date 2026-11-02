#!/usr/bin/env python3
"""Status, resources and headline scores for a comparator ROOT (results/comparators_<stamp>).

Usage: collect.py ROOT   -> writes ROOT/collected_<timestamp>.tsv and prints a progress table.
Resources come from sacct step 0 (the tool alone): AllocCPUS, Elapsed, MaxRSS. sacct samples RSS
every 30 s, so very short steps can show 0; qc.json max_rss (getrusage) is reported beside it.
Prototype held-out false-malignant rate is NA: it is 0 by construction (it never tests references).
Comparator 'split' recall comes from cn_state runs, not a statistical test, so it is labelled so.
"""
import datetime as dt
import glob
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


def sacct(job):
    out = subprocess.run(['sacct', '-j', job, '-P', '-n', '-o', 'JobID,State,AllocCPUS,Elapsed,MaxRSS'],
                         capture_output=True, text=True).stdout.splitlines()
    rows = {r.split('|')[0]: r.split('|') for r in out}
    top = rows.get(job, [job, 'UNKNOWN', '', '', ''])
    step = rows.get(f'{job}.0', ['', '', '', '', ''])
    return dict(state=top[1], cpus=step[2] or top[2], elapsed=step[3] or top[3], maxrss_sacct=step[4])


def main(root):
    root = Path(root)
    jobs = pd.read_csv(root / 'submitted_jobs.txt', sep=' ', header=None, names=['tool', 'dataset', 'job'], dtype=str)
    jobs = jobs.drop_duplicates(['tool', 'dataset'], keep='last')   # resubmissions replace failed jobs
    rows = []
    for r in jobs.itertuples():
        d = dict(tool=r.tool, dataset=r.dataset, job=r.job, **sacct(r.job))
        runs = sorted(glob.glob(str(root / r.dataset / f'{r.tool}_*' / 'qc.json')))
        if runs:
            q = json.loads(Path(runs[-1]).read_text())
            run = Path(runs[-1]).parent
            d['run'] = str(run)
            d['wall_seconds'] = q.get('wall_seconds')
            d['maxrss_gib_getrusage'] = q.get('max_rss_gib') or q.get('r_max_rss_gib')
            hr = q.get('holdout_malignant_rate')
            d['holdout_false_malignant'] = 'NA' if r.tool == 'prototype' else hr
            for sub in ['dna', 'eval']:
                s = run / sub / 'summary.json'
                if s.exists():
                    j = json.loads(s.read_text())
                    if sub == 'dna' and 'rho_log_ratio' in j:
                        d.update(rho_arms=j.get('rho_log_ratio'), altered_same_dir=j.get('altered_arms_same_direction'),
                                 altered=j.get('dna_altered_arms'), neutral_called=j.get('neutral_arms_called_altered'),
                                 neutral=j.get('neutral_arms'), tumour_clones=j.get('tumour_clones'),
                                 subclonal_split_cn_state_runs=j.get('subclonal_segments_split'),
                                 subclonal=j.get('subclonal_segments'),
                                 split_expected_by_chance=j.get('subclonal_split_fraction_expected_by_chance'))
                    elif sub == 'eval':
                        d.update(rho_arms_HCT116=(j.get('rho_arms_median') or {}).get('HCT116'),
                                 mixture_ari_line=j.get('mixture_ari_line_median'),
                                 mixture_clones=j.get('mixture_tumour_clones'),
                                 single_line_clones=j.get('single_line_tumour_clones'),
                                 single_line_ari_dose=j.get('single_line_ari_dose'),
                                 single_line_ari_batch=j.get('single_line_ari_batch'))
                    else:
                        d['dna_summary'] = json.dumps(j)[:300]
        rows.append(d)
    t = pd.DataFrame(rows)
    stamp = dt.datetime.now().strftime('%Y%m%d_%H%M%S')
    t.to_csv(root / f'collected_{stamp}.tsv', sep='\t', index=False)
    for tool, g in t.groupby('tool'):
        n = len(g)
        done = (g.state == 'COMPLETED').sum()
        run_ = (g.state == 'RUNNING').sum()
        pend = (g.state == 'PENDING').sum()
        bad = n - done - run_ - pend
        bar = '█' * round(20 * done / n) + '░' * (20 - round(20 * done / n))
        print(f'{tool:11s} {done:3d}/{n:<3d} {bar}  running {run_}  pending {pend}  failed/other {bad}')
    print(root / f'collected_{stamp}.tsv')


if __name__ == '__main__':
    main(sys.argv[1])
