"""Tahoe-100M (tahoebio/Tahoe-100M, CC0) extraction for claim 8 (docs/STATUS.md, "Claim 8 prespecified 2026-11-06").

Pinned to one dataset revision. Steps, all resumable (a finished step leaves its output and is skipped):
1. metadata: download obs/gene/cell_line/sample metadata; sha256 checked against the Hub LFS oid.
2. select: plate1 rows of obs_metadata with drug == 'DMSO_TF'. Writes the distinct pass_filter values (all rows,
   plate1, plate1 DMSO_TF) and the rule used, before any expression shard is touched. Rule: pass_filter == 'full';
   if that value does not exist the step stops (exit 3) for a decision.
3. footers: every expression shard's parquet footer by HTTP range read: rows and per-row-group min/max of plate
   and sample. Candidate shards: a row group whose plate range contains plate1 and whose sample range overlaps
   the DMSO_TF samples (or has no statistics).
4. shards: download candidates (sha256 vs LFS oid), keep plate1 DMSO_TF rows in plate1_dmso_tf.parquet. Shards
   with no hit are deleted. Matching is by content (plate, sample, BARCODE_SUB_LIB_ID), not by row order.
5. verify: every selected obs cell is found exactly once; counts written to verify.json. Exit 4 otherwise.
No expression value is read or summarised here; expression lists are carried through unchanged.
Usage: python -I fetch_tahoe.py OUTDIR
"""
import hashlib
import io
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import requests

REPO = 'tahoebio/Tahoe-100M'
REV = '349e2a55f5538a3e73b32a06c6bdd8e9e87011e8'
BASE = f'https://huggingface.co/datasets/{REPO}/resolve/{REV}/'
N_SHARDS = 3388
PLATE = 'plate1'
DRUG = 'DMSO_TF'
PASS = 'full'
META = ['obs_metadata', 'gene_metadata', 'cell_line_metadata', 'sample_metadata']
KEY = ['plate', 'sample', 'BARCODE_SUB_LIB_ID']


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def retry(f, *a, tries=6):
    for k in range(tries):
        try:
            return f(*a)
        except (requests.RequestException, OSError, AssertionError) as e:
            if k == tries - 1:
                raise
            log('retry', k + 1, repr(e)[:200])
            time.sleep(5 * 2 ** k)


def paths_info(paths):
    out = {}
    for i in range(0, len(paths), 200):
        r = requests.post(f'https://huggingface.co/api/datasets/{REPO}/paths-info/{REV}',
                          data={'paths': paths[i:i + 200]}, timeout=120)
        r.raise_for_status()
        for x in r.json():
            out[x['path']] = dict(size=x['size'], sha256=(x.get('lfs') or {}).get('oid'))
    return out


def download(path, dest, info):
    if dest.exists() and dest.stat().st_size == info['size']:
        h = sha256(dest)
        if h == info['sha256']:
            return h
    tmp = dest.with_suffix(dest.suffix + '.part')
    h = hashlib.sha256()
    with requests.get(BASE + path, stream=True, timeout=300) as r:
        r.raise_for_status()
        with open(tmp, 'wb') as f:
            for chunk in r.iter_content(1 << 22):
                f.write(chunk)
                h.update(chunk)
    got = h.hexdigest()
    assert tmp.stat().st_size == info['size'], f'{path}: size {tmp.stat().st_size} != {info["size"]}'
    assert got == info['sha256'], f'{path}: sha256 {got} != {info["sha256"]}'
    tmp.rename(dest)
    return got


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 24), b''):
            h.update(chunk)
    return h.hexdigest()


class RangeFile(io.RawIOBase):
    """Read-only seekable file over HTTP Range requests (enough for a parquet footer)."""

    def __init__(self, url, session):
        self.s, self.pos = session, 0
        r = session.head(url, allow_redirects=True, timeout=120)
        r.raise_for_status()
        self.url, self.size = r.url, int(r.headers['Content-Length'])

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off)
        return self.pos

    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0:
            return 0
        r = self.s.get(self.url, headers={'Range': f'bytes={self.pos}-{self.pos + n - 1}'}, timeout=300)
        r.raise_for_status()
        assert r.status_code == 206 and len(r.content) == n, (r.status_code, len(r.content), n)
        b[:n] = r.content
        self.pos += n
        return n


def step_metadata(out, manifest):
    d = out / 'metadata'
    d.mkdir(exist_ok=True)
    info = paths_info([f'metadata/{m}.parquet' for m in META])
    for m in META:
        p = f'metadata/{m}.parquet'
        h = retry(download, p, d / f'{m}.parquet', info[p])
        manifest.append(dict(path=p, local=str(d / f'{m}.parquet'), bytes=info[p]['size'], sha256=h, revision=REV))
        log('metadata', m, info[p]['size'], h)


def step_select(out):
    sel_path = out / 'selected_cells.tsv.gz'
    if sel_path.exists():
        return pd.read_csv(sel_path, sep='\t', dtype=str)
    samples = pd.read_parquet(out / 'metadata/sample_metadata.parquet')
    dmso = sorted(samples.loc[samples.plate.eq(PLATE) & samples.drug.eq(DRUG), 'sample'])
    log('DMSO_TF samples of', PLATE, dmso)
    pf = pq.ParquetFile(out / 'metadata/obs_metadata.parquet')
    cols = ['plate', 'sample', 'BARCODE_SUB_LIB_ID', 'BARCODE', 'sublibrary', 'drug', 'cell_line', 'cell_name',
            'pass_filter', 'tscp_count', 'gene_count', '__index_level_0__']
    counts, parts, offset = [], [], 0
    for g in range(pf.metadata.num_row_groups):
        t = pf.read_row_group(g, columns=['plate', 'drug', 'pass_filter'])
        t = pa.table({c: pc.cast(t[c], pa.string()) for c in t.column_names})
        n = t.num_rows
        vc = t.group_by(['plate', 'drug', 'pass_filter']).aggregate([([], 'count_all')]).to_pandas()
        counts.append(vc)
        m = pc.equal(t['plate'], PLATE)
        if pc.any(m).as_py():
            full = pf.read_row_group(g, columns=cols)
            keep = pc.and_(pc.equal(full['plate'], PLATE), pc.equal(full['drug'], DRUG))
            idx = pa.array(range(offset, offset + n))
            sub = full.append_column('obs_row', idx).filter(keep)
            parts.append(sub.to_pandas())
        offset += n
    vc = pd.concat(counts).groupby(['plate', 'drug', 'pass_filter'], observed=True)['count_all'].sum().reset_index()
    rows = [dict(scope='all', pass_filter=k, n=int(v)) for k, v in vc.groupby('pass_filter', observed=True).count_all.sum().items()]
    p1 = vc[vc.plate.eq(PLATE)]
    rows += [dict(scope=PLATE, pass_filter=k, n=int(v)) for k, v in p1.groupby('pass_filter', observed=True).count_all.sum().items()]
    p1d = p1[p1.drug.eq(DRUG)]
    rows += [dict(scope=f'{PLATE}_{DRUG}', pass_filter=k, n=int(v)) for k, v in p1d.groupby('pass_filter', observed=True).count_all.sum().items()]
    pd.DataFrame(rows).to_csv(out / 'pass_filter_values.tsv', sep='\t', index=False)
    log('pass_filter values\n' + pd.DataFrame(rows).to_string(index=False))
    log('obs rows', offset)
    sel = pd.concat(parts, ignore_index=True)
    for c in sel.columns:
        if isinstance(sel[c].dtype, pd.CategoricalDtype):
            sel[c] = sel[c].astype(str)
    if PASS not in set(sel.pass_filter):
        (out / 'selection_rule.txt').write_text(f'pass_filter value {PASS!r} absent in {PLATE} {DRUG}; stopped for a decision\n')
        log('STOP: no', PASS)
        sys.exit(3)
    assert set(sel['sample']) == set(dmso), (set(sel['sample']), dmso)
    assert not sel.duplicated(KEY).any(), 'key (plate, sample, BARCODE_SUB_LIB_ID) not unique'
    sel['selected'] = sel.pass_filter.eq(PASS)
    (out / 'selection_rule.txt').write_text(
        f'revision {REV}\nplate == {PLATE!r} and drug == {DRUG!r} (samples {dmso}); pass_filter == {PASS!r} '
        f'(values in this subset: {sorted(set(sel.pass_filter))}). Written {time.strftime("%Y-%m-%d %H:%M:%S")}, '
        'before any expression shard was read.\n')
    sel.to_csv(sel_path, sep='\t', index=False)
    log('selected', int(sel.selected.sum()), 'of', len(sel), 'plate1 DMSO_TF obs rows')
    return pd.read_csv(sel_path, sep='\t', dtype=str)


def footer(i, session):
    path = f'data/train-{i:05d}-of-{N_SHARDS:05d}.parquet'
    md = pq.ParquetFile(RangeFile(BASE + path, session)).metadata
    names = [md.schema.column(j).path for j in range(md.num_columns)]
    rows = []
    for g in range(md.num_row_groups):
        rg = md.row_group(g)
        r = dict(shard=i, path=path, row_group=g, rows=rg.num_rows)
        for c in ('plate', 'sample'):
            st = rg.column(names.index(c)).statistics
            ok = st is not None and st.has_min_max
            r[f'{c}_min'], r[f'{c}_max'] = (st.min, st.max) if ok else (None, None)
        rows.append(r)
    return rows


def step_footers(out):
    p = out / 'footer_index.tsv.gz'
    if p.exists():
        return pd.read_csv(p, sep='\t')
    session = requests.Session()
    with ThreadPoolExecutor(8) as ex:
        res = list(ex.map(lambda i: retry(footer, i, session), range(N_SHARDS)))
    d = pd.DataFrame([r for rows in res for r in rows])
    d.to_csv(p, sep='\t', index=False)
    log('footers', d.shard.nunique(), 'shards,', int(d.rows.sum()), 'expression rows')
    return d


def step_shards(out, foot, sel, manifest):
    dest = out / 'plate1_dmso_tf.parquet'
    hits_path = out / 'shard_hits.tsv'
    samples = sorted(set(sel['sample']))
    lo, hi = samples[0], samples[-1]
    na = foot[['plate_min', 'plate_max', 'sample_min', 'sample_max']].isna().any(axis=1)
    f = foot.fillna('')
    cand = foot[na | ((f.plate_min <= PLATE) & (f.plate_max >= PLATE) & (f.sample_min <= hi) & (f.sample_max >= lo))]
    shards = sorted(cand.path.unique())
    log('candidate shards', len(shards), 'of', foot.shard.nunique())
    if dest.exists() and hits_path.exists():
        return
    if len(shards) > 600:
        raise SystemExit(f'{len(shards)} candidate shards: row groups are not localised; stream-and-filter needed, stopping')
    info = paths_info(shards)
    sd = out / 'shards'
    sd.mkdir(exist_ok=True)
    writer, hits = None, []
    keyset = pa.array(samples)
    for path in shards:
        local = sd / Path(path).name
        h = retry(download, path, local, info[path])
        pf = pq.ParquetFile(local)
        t = pf.read(columns=['plate', 'sample'])
        m = pc.and_(pc.equal(t['plate'], PLATE), pc.is_in(t['sample'], value_set=keyset))
        n = int(pc.sum(m).as_py() or 0)
        hits.append(dict(path=path, bytes=info[path]['size'], sha256=h, rows=pf.metadata.num_rows, hits=n))
        if n:
            full = pf.read()
            rowidx = pa.array(range(full.num_rows), pa.int64())
            sub = full.append_column('shard', pa.array([path] * full.num_rows)).append_column('shard_row', rowidx).filter(m)
            if writer is None:
                writer = pq.ParquetWriter(dest.with_suffix('.parquet.part'), sub.schema, compression='zstd')
            writer.write_table(sub)
            manifest.append(dict(path=path, local=str(local), bytes=info[path]['size'], sha256=h, revision=REV))
        else:
            local.unlink()
        log(path, 'hits', n)
    writer.close()
    dest.with_suffix('.parquet.part').rename(dest)
    pd.DataFrame(hits).to_csv(hits_path, sep='\t', index=False)


def step_verify(out, sel):
    ext = pq.read_table(out / 'plate1_dmso_tf.parquet', columns=KEY + ['drug', 'cell_line_id']).to_pandas()
    e = ext.merge(sel[KEY + ['drug', 'cell_line', 'selected']], on=KEY, how='outer', indicator=True, suffixes=('_x', '_obs'))
    both = e._merge.eq('both')
    v = dict(
        obs_plate1_dmso_rows=len(sel), obs_selected=int(sel.selected.eq('True').sum()),
        extracted_rows=len(ext), extracted_key_duplicates=int(ext.duplicated(KEY).sum()),
        matched=int(both.sum()), only_in_extracted=int(e._merge.eq('left_only').sum()),
        only_in_obs=int(e._merge.eq('right_only').sum()),
        selected_missing=int((e._merge.eq('right_only') & e.selected.eq('True')).sum()),
        drug_mismatch=int((both & e.drug_x.ne(e.drug_obs)).sum()),
        cell_line_mismatch=int((both & e.cell_line_id.ne(e.cell_line)).sum()),
    )
    (out / 'verify.json').write_text(json.dumps(v, indent=1))
    log('verify', v)
    bad = v['extracted_key_duplicates'] or v['selected_missing'] or v['drug_mismatch'] or v['cell_line_mismatch']
    if bad:
        sys.exit(4)


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    manifest = []
    step_metadata(out, manifest)
    sel = step_select(out)
    foot = step_footers(out)
    step_shards(out, foot, sel, manifest)
    old = pd.read_csv(out / 'manifest.tsv', sep='\t') if (out / 'manifest.tsv').exists() else pd.DataFrame()
    pd.concat([old, pd.DataFrame(manifest)]).drop_duplicates('path', keep='last').to_csv(out / 'manifest.tsv', sep='\t', index=False)
    step_verify(out, sel)
    log('done')


if __name__ == '__main__':
    main(*sys.argv[1:])
