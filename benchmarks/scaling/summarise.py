"""Summarise the compute scaling benchmark (run_one.sbatch outputs under ROOT/runs/<tool>_<N>_rep<k>/).

Per run: wall seconds from /usr/bin/time -v (whole pipeline: input reading, the tool, our post-processing for
comparators), peak memory as sacct MaxRSS of the batch step (the whole process tree; /usr/bin/time reports one
process, which undercounts R tools that fork workers), CPU utilisation, the tool's own core seconds where its qc.json records them
(comparators: r_seconds or inference_seconds), node, and status from sacct (COMPLETED, TIMEOUT, OUT_OF_MEMORY,
FAILED). Per tool and size: median over repeats. All runs: cpu partition, --constraint=v6, 16 CPUs.
Usage: summarise.py ROOT
"""
import json
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd


def parse_time(path):
    out = {}
    if not path.exists():
        return out
    for line in path.read_text(errors='replace').splitlines():
        line = line.strip()
        if line.startswith('Elapsed (wall clock) time'):
            v = line.rsplit(': ', 1)[1].split(':')
            out['wall_s'] = sum(float(x) * 60 ** i for i, x in enumerate(reversed(v)))
        elif line.startswith('Maximum resident set size'):
            out['peak_rss_gib'] = int(line.rsplit(': ', 1)[1]) / 1024 ** 2
        elif line.startswith('Percent of CPU this job got'):
            out['cpu_percent'] = float(line.rsplit(': ', 1)[1].rstrip('%'))
    return out


def main(root):
    root = Path(root)
    rows = []
    for d in sorted((root / 'runs').glob('*_rep*')):
        m = re.fullmatch(r'(.+)_(\d+)_rep(\d+)', d.name)
        tool, n, rep = m.group(1), int(m.group(2)), int(m.group(3))
        job = dict(l.split('\t', 1) for l in (d / 'job.tsv').read_text().splitlines() if '\t' in l)
        state = subprocess.run(['sacct', '-j', job['slurm_job'], '-X', '-n', '-o', 'State%20'],
                               capture_output=True, text=True).stdout.strip()
        rss = subprocess.run(['sacct', '-j', f"{job['slurm_job']}.batch", '-n', '-o', 'MaxRSS%30', '--noconvert'],
                             capture_output=True, text=True).stdout.strip()
        r = dict(tool=tool, cells=n, rep=rep, host=job.get('host'), exit=job.get('exit'), state=state)
        r.update(parse_time(d / 'time.txt'))
        r['time_rss_gib'] = r.pop('peak_rss_gib', None)
        rss = rss.split()[0] if rss.split() else ''
        sacct_gib = int(rss) / 1024 ** 3 if rss.isdigit() else None
        # sacct samples every 30 s (misses short peaks); /usr/bin/time sees one process. Take the larger.
        r['peak_rss_gib'] = max(x for x in (sacct_gib, r['time_rss_gib'], 0.0) if x is not None)
        qc = list(d.glob('*/qc.json'))
        if qc:
            q = json.loads(qc[0].read_text())
            r['core_s'] = q.get('r_seconds', q.get('inference_seconds'))
        rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(root / 'runs.tsv', sep='\t', index=False)
    ok = t[t.state.eq('COMPLETED') & t.exit.eq('0')]
    s = ok.groupby(['tool', 'cells']).agg(repeats=('rep', 'size'), wall_min=('wall_s', lambda x: x.median() / 60),
                                          peak_rss_gib=('peak_rss_gib', 'median'),
                                          cpu_percent=('cpu_percent', 'median')).round(2).reset_index()
    bad = t[~t.index.isin(ok.index)][['tool', 'cells', 'rep', 'state', 'exit']]
    s.to_csv(root / 'summary.tsv', sep='\t', index=False)
    bad.to_csv(root / 'not_finished.tsv', sep='\t', index=False)
    print(s.to_string(index=False))
    if len(bad):
        print('\nnot finished:\n' + bad.to_string(index=False))


if __name__ == '__main__':
    main(*sys.argv[1:])
