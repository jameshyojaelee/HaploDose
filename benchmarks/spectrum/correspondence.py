"""Metadata-only DNA/RNA correspondence; anatomical names never prove same biopsy.

Usage: correspondence.py RNA_SAMPLES_TSV DNA_CELL_INFO_DIR OUT
Reads only the supplied RNA sample manifest and small *_cell_info.csv.gz tables.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


def site_name(value):
    return re.sub(r'[^a-z0-9]', '', str(value).lower()) if pd.notna(value) else ''


def rna_site(sample):
    parts = re.split(r'_(?:CD45N|CD45P|UNSORTED)_', str(sample), maxsplit=1)
    return parts[1] if len(parts) == 2 else ''


def correspondence(samples, dna_dir):
    if samples.sample_title.isna().any() or samples.sample_title.duplicated().any():
        raise ValueError('RNA sample identifiers must be nonmissing and unique')
    rows = []
    for patient, group in samples.groupby('patient', sort=True):
        path = Path(dna_dir) / f'{patient}_cell_info.csv.gz'
        ci = pd.read_csv(path, dtype={'cell_id': str, 'leaf_id': str}) if path.exists() else pd.DataFrame()
        if len(ci) and (ci.cell_id.isna().any() or ci.cell_id.duplicated().any()):
            raise ValueError(f'duplicated or missing DNA cell identifiers: {patient}')
        sites = ci.tumor_site.fillna('').unique() if len(ci) else ['']
        for sample in group.sample_title:
            rs = rna_site(sample)
            for ds in sites:
                subset = ci[ci.tumor_site.fillna('').eq(ds)] if len(ci) else ci
                n = len(subset)
                tier = ('missing_site' if not site_name(rs) or not site_name(ds) else
                        'same_anatomical_label_only' if site_name(rs) == site_name(ds) else 'unresolved_site_mapping')
                rows.append(dict(patient=patient, rna_specimen=sample, rna_site=rs, dna_site=ds,
                                 match_tier=tier, same_biopsy_verified=False, time_match='unknown',
                                 evidence=str(path) if path.exists() else '', dna_cells=n,
                                 dna_published_clones=int(subset.leaf_id.nunique()) if n else 0,
                                 unseen_clone_prevalence_upper95=float(-np.expm1(np.log(0.05) / n)) if n else np.nan,
                                 prevalence_bound_assumption='independent representative DNA-cell sampling'))
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('samples', type=Path)
    ap.add_argument('dna_dir', type=Path)
    ap.add_argument('out', type=Path)
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=False)
    result = correspondence(pd.read_csv(a.samples, sep='\t'), a.dna_dir)
    result.to_csv(a.out / 'correspondence.tsv', sep='\t', index=False)
    (a.out / 'summary.json').write_text(json.dumps(dict(
        schema_version=2, rows=len(result), same_biopsy_verified=0,
        tiers=result.match_tier.value_counts().to_dict() if len(result) else {},
        limit='Site-label agreement is not same-biopsy provenance. DNA/RNA cells are separate; '
              'sampling bounds do not account for capture, QC or anatomical bias.'), indent=2) + '\n')


if __name__ == '__main__':
    main()
