"""Haplodose bulk inputs for the REBC-THYR benchmark (prespecified claim 7, docs/STATUS.md "Claim 7 prespecified
2026-11-06"). benchmarks/cptac_bulk/make_input.py unchanged except the cohort map: one group THYR.
Amendment 1 (2026-11-06, before any claim 7 score was computed): REBC RNA files carry sequencing-centre sample IDs
(SC*_merged), never the ASCAT tumour sample ID, so the claim 6 sample match finds 5 of 383 cases. Cases with no
sample-ID match fall back to a case-level match: the case's Primary Tumor STAR - Counts files, several -> smallest
file_id. Recorded per case in selection.tsv (column match = sample | case).
Usage: make_input.py DATA OUTDIR (outputs as the CPTAC builder).
"""
import importlib.util
import sys

import numpy as np
import pandas as pd
from pathlib import Path

H = Path('/gpfs/commons/home/jameslee/Haplodose')
_spec = importlib.util.spec_from_file_location('cptac_make_input', H / 'benchmarks/cptac_bulk/make_input.py')
cmi = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cmi)


def cohort(site, disease):
    """REBC-THYR group from GDC case primary_site and disease_type; None = not in claim 7."""
    if (site or '').lower() == 'thyroid gland' and (disease or '').lower() == 'adenomas and adenocarcinomas':
        return 'THYR'
    return None


_select_sample = cmi.select


def select(meta):
    """Claim 6 sample-ID match, then amendment 1 case-level fallback for unmatched cases."""
    tum, nor, qc = _select_sample(meta)
    qc['match'] = np.where(qc.included, 'sample', '')
    m = meta.copy()
    m['cohort'] = [cohort(s, d) for s, d in zip(m.primary_site, m.disease_type)]
    m = m[m.cohort.notna()]
    asc = m[m.kind.eq('ascat')].sort_values('file_id').drop_duplicates('case').set_index('case')
    rna = m[m.kind.eq('rna')]
    pt = rna[rna.sample_types.str.split(';').map(lambda x: all(v == 'Primary Tumor' for v in x))].sort_values('file_id')
    rows = []
    for case in qc.case[~qc.included]:
        c = pt[pt.case.eq(case)]
        if len(c):
            r, a = c.iloc[0], asc.loc[case]
            rows.append(dict(case=case, cohort=a.cohort, rna_file_id=r.file_id, rna_path=r.path,
                             pooled=len(r.sample_ids.split(';')) > 1, ascat_file_id=a.file_id, ascat_path=a.path))
            qc.loc[qc.case.eq(case), ['n_rna_candidates', 'included', 'match']] = [len(c), True, 'case']
    tum = pd.concat([tum, pd.DataFrame(rows)], ignore_index=True)
    return tum, nor, qc


cmi.cohort = cohort
cmi.COHORTS = ['THYR']
cmi.select = select

if __name__ == '__main__':
    cmi.main(sys.argv[1], sys.argv[2])
