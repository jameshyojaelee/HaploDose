"""Locked three-line RNA/DNA development diagnostic, run with the installed wheel.

No same-cell membership or absolute DNA-copy claim. Raw selected counts are checked
against GEO, annotation against a supplied hg38 GTF, before any model is fitted.
Example resources: cpu, 8 CPUs, 32G, 90 hours (single job; not an array).
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import gzip
import hashlib
import io
from importlib import metadata
import inspect
import json
from pathlib import Path
import re
import shutil
import tarfile

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from sklearn.metrics import roc_auc_score

LINES = ('SNU-601', 'MKN-45', 'NUGC-4')


class IdentityGap(ValueError):
    def __init__(self, gaps):
        self.gaps = gaps
        super().__init__('ambiguous retained raw symbols; source rows are not bound to unique GTF loci')


def record(path):
    path = Path(path)
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            h.update(block)
    return dict(path=str(path.resolve()), bytes=path.stat().st_size, sha256=h.hexdigest())


def read_input(directory):
    cells = pd.read_csv(directory / 'cells.tsv', sep='\t', keep_default_na=False)
    positions = pd.read_csv(directory / 'positions.tsv', sep='\t')
    with h5py.File(directory / 'counts.h5') as h:
        genes = h['genes'].asstr()[:]
        ids = h['cell_ids'].asstr()[:]
        X = sparse.csr_matrix((h['X/data'][:], h['X/indices'][:], h['X/indptr'][:]),
                              shape=(len(ids), len(genes)))
    if not np.array_equal(ids, cells.cell_id) or not np.array_equal(genes, positions.gene):
        raise ValueError('prepared matrix/metadata/gene order mismatch')
    if cells.cell_id.duplicated().any() or pd.Index(genes).duplicated().any():
        raise ValueError('duplicate prepared identities')
    if np.any(~np.isfinite(X.data) | (X.data < 0) | (X.data != np.floor(X.data))):
        raise ValueError('non-raw count values')
    return cells, positions, X


def annotation_gate(inputs, gtf):
    """Exact coordinates, no guessed conversion. GTF is one-based closed."""
    wanted = set().union(*(set(p.gene) for _, p, _ in inputs.values()))
    loci = {g: set() for g in wanted}
    identities = {g: {} for g in wanted}
    opener = gzip.open if str(gtf).endswith('.gz') else open
    assembly = False
    with opener(gtf, 'rt') as stream:
        for row in stream:
            if row.startswith('#'):
                assembly |= 'GRCh38' in row or 'hg38' in row
                continue
            fields = row.rstrip().split('\t')
            if len(fields) != 9 or fields[2] != 'gene':
                continue
            hit = re.search(r'gene_name "([^"]+)"', fields[8])
            if hit and hit[1] in wanted:
                chrom = fields[0] if fields[0].startswith('chr') else 'chr' + fields[0]
                locus = (chrom, int(fields[3]), int(fields[4]))
                loci[hit[1]].add(locus)
                identifier = re.search(r'gene_id "([^"]+)"', fields[8])
                if identifier:
                    identities[hit[1]].setdefault(locus, set()).add(identifier[1].split('.')[0])
    if not assembly:
        raise ValueError('GTF header does not independently declare hg38/GRCh38')
    for line, (_, pos, _) in inputs.items():
        bad = [r.gene for r in pos.itertuples() if
               (r.chromosome, int(r.start), int(r.end)) not in loci[r.gene]]
        if bad:
            raise ValueError(f'{line}: annotation mismatch for {len(bad)} genes; first={bad[0]}')
    return {r.gene: identities[r.gene].get((r.chromosome, int(r.start), int(r.end)), set())
            for _, pos, _ in inputs.values() for r in pos.itertuples()}


def raw_gate(inputs, tar_path, params, identities=None, diagnostics=None, identity_only=False,
             derive=False, canonical=False):
    """Countercheck every selected query/reference row, keeping original membership."""
    checked = {line: np.zeros(len(cells), bool) for line, (cells, _, _) in inputs.items()}
    gaps = []
    mappings = []
    deltas = []
    with tarfile.open(tar_path) as archive:
        for prefix, spec in sorted(params['libraries'].items()):
            library = prefix.split('_', 1)[1]
            selected = {line: np.flatnonzero(cells.library.eq(library))
                        for line, (cells, _, _) in inputs.items()}
            if not any(len(v) for v in selected.values()):
                continue
            def payload(suffix):
                return gzip.decompress(archive.extractfile(f'{prefix}_{suffix}').read())
            features = pd.read_csv(io.BytesIO(payload('genes.tsv.gz')), sep='\t', header=None)
            genes = features.iloc[:, -1].astype(str)
            retained = set().union(*(set(inputs[line][1].gene) for line, indices in selected.items() if len(indices)))
            ambiguous = sorted(set(genes[genes.duplicated(keep=False)]) & retained)
            unresolved = []
            canonical_rows = {}
            for symbol in ambiguous:
                raw = features.loc[genes.eq(symbol)].astype(str)
                ids = raw.iloc[:, 0].str.split('.').str[0].tolist() if raw.shape[1] > 1 else []
                target = sorted((identities or {}).get(symbol, set()))
                matches = [i for i, identifier in enumerate(ids) if identifier in target]
                resolved = len(target) == 1 and matches == [0]
                unique = len(target) == 1 and len(matches) == 1
                if unique:
                    canonical_rows[symbol] = int(raw.index[matches[0]])
                mappings.append(dict(library=library, symbol=symbol, raw_rows=raw.values.tolist(),
                                     canonical_locus_gtf_ids=target, matching_raw_row_offsets=matches,
                                     historical_first_row_certified=resolved,
                                     reason='unique canonical ID matches historical first row' if resolved else
                                     'requires new input derivative' if len(target) == 1 and len(matches) == 1 else
                                     'unmapped or ambiguous source ID/locus'))
                # The sole observed correction is locked before fitting and is
                # determined from RNA identities alone, never DNA outcomes.
                permitted_repair = (unique and symbol == 'PRICKLE4' and library in
                                    ['HGC-27_1', 'HGC-27_2', 'SNU-16_1', 'SNU-16_2', 'SNU-638'])
                if not resolved and not ((derive or canonical) and permitted_repair):
                    unresolved.append(symbol)
            if unresolved:
                gaps.append(dict(library=library, symbols=unresolved,
                                 raw_feature_rows=features.loc[genes.isin(unresolved)].astype(str).values.tolist(),
                                 limit='historical first row is not uniquely bound to canonical GTF locus/ID'))
                continue
            if identity_only:
                continue
            barcodes = pd.read_csv(io.BytesIO(payload('barcodes.tsv.gz')), sep='\t', header=None)[0].astype(str)
            ids = pd.Index(library + '_' + barcodes)
            if ids.duplicated().any():
                raise ValueError('duplicate GEO barcode identity')
            X = sparse.csr_matrix(mmread(io.BytesIO(payload('matrix.mtx.gz'))).T)
            if np.any(~np.isfinite(X.data) | (X.data < 0) | (X.data != np.floor(X.data))):
                raise ValueError('GEO source is not integer raw counts')
            first = ~genes.duplicated()
            for line, indices in selected.items():
                if not len(indices):
                    continue
                cells, pos, prepared = inputs[line]
                rows = ids.get_indexer(cells.iloc[indices].cell_id)
                columns = pd.Index(genes[first]).get_indexer(pos.gene)
                if np.any(rows < 0) or np.any(columns < 0):
                    raise ValueError('selected identity absent from GEO')
                columns = np.flatnonzero(first)[columns]
                original_columns = columns.copy()
                if canonical:
                    columns = np.array([canonical_rows.get(g, col) for g, col in zip(pos.gene, columns)])
                if (X[rows][:, columns] != prepared[indices]).nnz:
                    raise ValueError(f'{line}/{library}: selected counts differ from GEO')
                if derive:
                    for symbol, raw_row in canonical_rows.items():
                        column = int(np.flatnonzero(pos.gene.eq(symbol))[0])
                        if raw_row == original_columns[column]:
                            continue
                        old = prepared[indices, column].toarray().ravel().astype(np.int64)
                        new = X[rows, raw_row].toarray().ravel().astype(np.int64)
                        difference = new - old
                        changed = difference != 0
                        delta = sparse.csr_matrix((difference[changed],
                                                   (indices[changed], np.full(changed.sum(), column))),
                                                  shape=prepared.shape)
                        prepared = prepared + delta
                        for idx, before, after in zip(indices[changed], old[changed], new[changed]):
                            deltas.append(dict(line=line, library=library, gene=symbol,
                                               cell_id=cells.cell_id.iloc[idx], old=int(before), new=int(after)))
                    inputs[line] = (cells, pos, prepared)
                if np.any(np.asarray(X[rows].sum(axis=1)).ravel() < params['min_umi']):
                    raise ValueError('selected cell violates original raw-depth filter')
                part = cells.iloc[indices]
                if not part.patient.eq(spec['line']).all():
                    raise ValueError('patient/library mismatch')
                if spec['line'] == line:
                    valid = part.role.eq('query') & part.specimen.eq(library)
                else:
                    valid = (part.role.eq('reference') & part.ref_group.eq(spec['line']) &
                             part.ref_split.isin(['fit', 'holdout']) & part.specimen.eq('ref_' + spec['line']))
                if not valid.all():
                    raise ValueError('reference/query membership provenance mismatch')
                checked[line][indices] = True
    if diagnostics is not None:
        diagnostics['raw_duplicate_mapping'] = mappings
        if derive:
            diagnostics['count_deltas'] = deltas
    if gaps:
        raise IdentityGap(gaps)
    if identity_only:
        return {}
    if not all(v.all() for v in checked.values()):
        raise ValueError('some selected cells were not counterchecked')
    return {line: int(v.sum()) for line, v in checked.items()}


def supported_scores(rna, dna, line, dest):
    """Compare RNA bin spans only where corrected DNA has callable support.

    RNA bin spans are not base-resolution measurements. Arm coverage is exported;
    an arm with incomplete RNA clone coverage is missing, never a zero effect.
    """
    bins = pd.read_csv(rna / 'bins.tsv', sep='\t')
    profiles = pd.read_csv(rna / 'profiles.tsv.gz', sep='\t')
    cells = pd.read_csv(rna / 'cells.tsv.gz', sep='\t', dtype={'clone': str})
    db = pd.read_csv(dna / 'clone_bins.tsv', sep='\t')
    qc = json.loads((dna / 'qc.json').read_text())
    if qc.get('schema_version') != 2:
        raise ValueError('corrected schema-2 DNA truth required')
    weights = pd.Series(qc['clones'], dtype=float)
    dna_names = list(weights.index)
    output, rows = [], []
    query = cells[cells.role.eq('query') & cells.patient.eq(line)]
    for specimen, group in query.groupby('specimen'):
        sizes = group.clone.dropna().value_counts()
        sizes = sizes[~sizes.index.isin(['NA', 'unresolved', 'unassigned', 'excluded'])]
        sizes = sizes[sizes >= max(20, .02 * len(group))]
        selected = profiles[profiles.specimen.eq(specimen) & profiles.clone.isin(sizes.index)]
        pivot = selected.pivot(index='bin', columns='clone', values='log_ratio').reindex(columns=sizes.index)
        overlap_rows = []
        for k, d in db.iterrows():
            covered = bins[(bins.chromosome == d.chromosome) & (bins.arm == d.arm) &
                           (bins.start < d.end) & (bins.end > d.start)]
            cuts = sorted({int(x) for b in covered.itertuples() for x in
                           [max(b.start, d.start), min(b.end, d.end)]})
            for lo, hi in zip(cuts[:-1], cuts[1:]):
                active = covered.index[(covered.start < hi) & (covered.end > lo)]
                # Neighboring gene-bin spans can overlap. Share each elementary
                # interval equally across covering predictions; count its bp once.
                for j in active:
                    overlap_rows.append((j, k, (hi - lo) / len(active)))
        support = pd.DataFrame(overlap_rows, columns=['rna_bin', 'dna_bin', 'bp'])
        for arm, d in db.groupby('arm'):
            overlap = support[support.dna_bin.isin(d.index)]
            if overlap.empty or not len(sizes):
                continue
            rr = pivot.reindex(overlap.rna_bin).to_numpy(float)
            good = np.isfinite(rr).all(axis=1)
            overlap, rr = overlap.loc[good], rr[good]
            if overlap.empty:
                continue
            dd = db.loc[overlap.dna_bin, dna_names].to_numpy(float)
            w = overlap.bp.to_numpy(float)
            # Keep each pair fixed across an arm; changing pairs at every bin
            # would create arm truth unsupported by any individual clone pair.
            rd = max((abs(float(np.average(rr[:, i] - rr[:, j], weights=w)))
                      for i in range(rr.shape[1]) for j in range(i + 1, rr.shape[1])), default=0.)
            frac = max((float(np.average(np.abs(dd[:, i] - dd[:, j]) >= 1, weights=w))
                        for i in range(dd.shape[1]) for j in range(i + 1, dd.shape[1])), default=0.)
            r = dict(line=line, specimen=specimen, arm=arm, rna_groups=len(sizes),
                     dna_groups=len(dna_names), dna_callable_bp=int((d.end-d.start).sum()),
                     overlapping_span_bp=int(round(w.sum())), rna_max_pair_delta=rd,
                     dna_max_fixed_pair_fraction_ge1copy=frac, dna_variable_arm=frac >= .3)
            output.append(r)
        frame = pd.DataFrame([r for r in output if r['specimen'] == specimen])
        defined = len(frame) and frame.dna_variable_arm.nunique() == 2 and len(sizes) > 1
        rows.append(dict(line=line, specimen=specimen, query_cells=len(group), assigned_cells=int(sizes.sum()),
                         rna_groups=len(sizes), dna_groups=len(dna_names), supported_arms=len(frame),
                         profile_difference_auroc=float(roc_auc_score(frame.dna_variable_arm, frame.rna_max_pair_delta))
                         if defined else None,
                         limitation='arm-level unsigned profile differences; no clone correspondence or membership truth'))
    pd.DataFrame(output).to_csv(dest / 'supported_arms.tsv', sep='\t', index=False)
    pd.DataFrame(rows).to_csv(dest / 'scores.tsv', sep='\t', index=False)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-root', type=Path, required=True)
    ap.add_argument('--truth-root', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--annotation-gtf', type=Path)
    ap.add_argument('--identity-only', action='store_true', help='check source identity/annotation only; do not fit or certify counts')
    ap.add_argument('--derive-inputs', action='store_true', help='write only the independently ID-established PRICKLE4 correction to new inputs')
    ap.add_argument('--input-root', type=Path, help='new input derivative with a bound derivative.json certificate')
    a = ap.parse_args()
    root, truth, out = a.data_root.resolve(), a.truth_root.resolve(), a.out.resolve()
    if out.exists() or any(p == out or p in out.parents for p in [root / 'data', root / 'results', truth]):
        raise ValueError('output must be new and outside immutable data/results/truth')
    out.mkdir(parents=True)
    status = dict(schema_version=2, status='blocked', biological_n=3, role='locked development diagnostic',
                  limits=['separate RNA/DNA cells; no matched memberships',
                          'other cancer lines are references, not normal DNA; cell state and reference confounding remain',
                          'low-depth single DNA clusters do not establish clonality',
                          'RNA span overlap is not base-resolution measurement'])
    try:
        import haplodose
        from haplodose.clones import CloneParams
        from haplodose.pipeline import run
        module_path = Path(haplodose.__file__).resolve()
        if (root / 'python') in module_path.parents:
            raise ValueError('import resolves to live repository source, not installed wheel')
        if a.identity_only and a.derive_inputs:
            raise ValueError('identity-only and derive-inputs are separate modes')
        if a.derive_inputs and a.input_root:
            raise ValueError('derive only from original locked inputs')
        source = a.input_root.resolve() if a.input_root else root / 'results/andor_20261101_1224/rna_input_v2'
        if a.input_root:
            certificate = json.loads((source / 'derivative.json').read_text())
            if certificate.get('status') != 'raw_verified_derivative' or certificate.get('lines') != list(LINES):
                raise ValueError('invalid derivative certificate')
            for name, expected in certificate['output_sha256'].items():
                if record(source / name)['sha256'] != expected:
                    raise ValueError(f'derivative certificate mismatch: {name}')
        params = json.loads((source / 'params.json').read_text())
        inputs = {line: read_input(source / line) for line in LINES}
        files = [source / 'params.json', Path(__file__),
                 Path(__file__).resolve().parents[1] / 'andor/prepare_rna.py',
                 module_path, Path(inspect.getfile(CloneParams)), Path(inspect.getfile(run))]
        files += [source / line / name for line in LINES for name in ['counts.h5', 'cells.tsv', 'positions.tsv']]
        tar_path = root / 'data/external/andor2020_gse142750/GSE142750_RAW.tar'
        files.append(tar_path)
        status['files'] = [record(p) for p in files]
        if a.annotation_gtf is None:
            raise ValueError('--annotation-gtf required to establish assembly and exact one-based coordinates')
        status['annotation_file'] = record(a.annotation_gtf)
        identities = annotation_gate(inputs, a.annotation_gtf)
        files.append(a.annotation_gtf)
        status['annotation_provenance'] = ('Exact coordinate/content match to supplied hg38 GTF; '
                                           'does not recover historical execution provenance. '
                                           'Historical builder GRCh38-2020-A label is unverified.')
        status['raw_counterchecked_cells'] = raw_gate(inputs, tar_path, params, identities, status, a.identity_only,
                                                     derive=a.derive_inputs, canonical=bool(a.input_root))
        if a.derive_inputs:
            derivative = out / 'input'
            derivative.mkdir()
            original_inputs = {line: read_input(source / line) for line in LINES}
            for line, (cells, pos, X) in inputs.items():
                original = original_inputs[line][2]
                difference = (X - original).tocoo()
                if any(pos.gene.iloc[j] != 'PRICKLE4' or cells.role.iloc[i] != 'reference' or
                       cells.library.iloc[i] not in ['HGC-27_1', 'HGC-27_2', 'SNU-16_1', 'SNU-16_2', 'SNU-638']
                       for i, j in zip(difference.row, difference.col)):
                    raise ValueError('derivative changed a nonaffected count')
                if np.any(X.data < 0) or np.any(X.data > np.iinfo(np.int32).max):
                    raise ValueError('invalid corrected count')
                dest = derivative / line
                dest.mkdir()
                shutil.copy2(source / line / 'counts.h5', dest / 'counts.h5')
                with h5py.File(dest / 'counts.h5', 'r+') as h:
                    for key, values in [('data', X.data.astype(np.int32)), ('indices', X.indices.astype(np.int32)),
                                        ('indptr', X.indptr.astype(np.int64))]:
                        del h['X/' + key]
                        h.create_dataset('X/' + key, data=values, compression='gzip' if key != 'indptr' else None)
                for name in ['cells.tsv', 'positions.tsv']:
                    shutil.copy2(source / line / name, dest / name)
                    if record(dest / name)['sha256'] != record(source / line / name)['sha256']:
                        raise ValueError('metadata changed during derivative creation')
            shutil.copy2(source / 'params.json', derivative / 'params.json')
            reread = {line: read_input(derivative / line) for line in LINES}
            raw_gate(reread, tar_path, params, identities, canonical=True)
            status.update(status='raw_verified_derivative', lines=list(LINES),
                          output_sha256={str(p.relative_to(derivative)): record(p)['sha256']
                                         for p in sorted(derivative.glob('*/*')) if p.is_file()},
                          correction='PRICKLE4 reference entries uniquely mapped to ENSG00000278224; no other counts changed')
            status['output_sha256']['params.json'] = record(derivative / 'params.json')['sha256']
            (derivative / 'derivative.json').write_text(json.dumps(status, indent=2) + '\n')
            return
        if a.identity_only:
            status['status'] = 'identity_checked'
            status['counts_verified'] = False
            return
        settings = CloneParams(split_null_draws=0)
        status.update(package_version=metadata.version('haplodose'), installed_module=str(module_path),
                      clone_params=asdict(settings), bin_genes=10, chunk=4096, exclude_ribosomal=True,
                      counts_contract=dict(counts_kind='raw', genome_build='hg38',
                                           coordinate_system='one-based-closed', gene_namespace='symbol'),
                      configuration_change='selected-parent null disabled because current gene-aware null unsupported; '
                      'not a pure accuracy regression against historical null-filtered runs')
        status['files'] = [record(p) for p in files]
        (out / 'provenance.json').write_text(json.dumps(status, indent=2) + '\n')
        scores = []
        for line in LINES:
            directory, dest = source / line, out / line
            run(directory / 'counts.h5', directory / 'cells.tsv', directory / 'positions.tsv', dest,
                params=settings, bin_genes=10, chunk=4096, exclude_ribosomal=True,
                counts_kind='raw', genome_build='hg38', coordinate_system='one-based-closed', gene_namespace='symbol')
            dna = truth / ('andor_full' if line == 'SNU-601' else f'andor_low_{line}') / line
            status.setdefault('truth_files', []).extend(record(dna / f) for f in ['qc.json', 'clone_bins.tsv'])
            scores.extend(supported_scores(dest, dna, line, dest))
        pd.DataFrame(scores).to_csv(out / 'scores.tsv', sep='\t', index=False)
        status['status'] = 'completed'
    except Exception as exc:
        status['error'] = f'{type(exc).__name__}: {exc}'
        if isinstance(exc, IdentityGap):
            status['raw_identity_gaps'] = exc.gaps
        raise
    finally:
        (out / 'summary.json').write_text(json.dumps(status, indent=2) + '\n')


if __name__ == '__main__':
    main()
