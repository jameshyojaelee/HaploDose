"""GDC metadata for the REBC-THYR claim 7 files (no data content read). Same output and filters as
benchmarks/cptac_bulk/gdc_manifest.py with project REBC-THYR and the claim 7 cohort map.
Usage: gdc_manifest.py OUTDIR
"""
import importlib.util
import json
import sys
import urllib.request
from pathlib import Path

import pandas as pd

H = Path('/gpfs/commons/home/jameslee/Haplodose')


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gm = _load('cptac_gdc_manifest', H / 'benchmarks/cptac_bulk/gdc_manifest.py')
mi = _load('rebc_make_input', H / 'benchmarks/rebc_bulk/make_input.py')


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    proj = gm.f_in('cases.project.project_id', ['REBC-THYR'])
    rna = gm.files({'op': 'and', 'content': [proj, gm.f_in('data_type', ['Gene Expression Quantification']),
                                             gm.f_in('analysis.workflow_type', ['STAR - Counts']),
                                             gm.f_in('cases.samples.sample_type', ['Primary Tumor', 'Solid Tissue Normal'])]})
    asc = gm.files({'op': 'and', 'content': [proj, gm.f_in('data_type', ['Allele-specific Copy Number Segment']),
                                             gm.f_in('analysis.workflow_type', ['AscatNGS'])]})
    t = pd.DataFrame(gm.rows(rna, 'rna') + gm.rows(asc, 'ascat'))
    t['cohort'] = [mi.cohort(s, d) for s, d in zip(t.primary_site, t.disease_type)]
    t = t[t.cohort.notna()].sort_values(['kind', 'file_id'])
    t.to_csv(out / 'files.tsv', sep='\t', index=False)
    st = json.load(urllib.request.urlopen('https://api.gdc.cancer.gov/status'))
    (out / 'gdc_status.json').write_text(json.dumps(st, indent=1))
    print(t.groupby(['cohort', 'kind']).size().unstack(), '\nGB', t.file_size.astype(float).sum() / 1e9,
          '\naccess', t.access.unique())


if __name__ == '__main__':
    main(sys.argv[1])
