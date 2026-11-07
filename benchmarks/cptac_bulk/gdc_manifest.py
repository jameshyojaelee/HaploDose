"""GDC metadata for the CPTAC-3 claim 6 files (no data content read). Writes OUTDIR/files.tsv in the format
make_input.py expects, with the GDC md5sum and file_size used by fetch.sbatch.
  rna    CPTAC-3 Gene Expression Quantification, workflow STAR - Counts, any sample Primary Tumor or Solid Tissue Normal
  ascat  CPTAC-3 Allele-specific Copy Number Segment, workflow AscatNGS
Only files whose case maps to one of the 9 claim 6 cohorts (make_input.cohort) are listed.
Usage: gdc_manifest.py OUTDIR
"""
import importlib.util
import json
import sys
import urllib.request
from pathlib import Path

import pandas as pd

H = Path('/gpfs/commons/home/jameslee/Haplodose')
_spec = importlib.util.spec_from_file_location('cptac_make_input', H / 'benchmarks/cptac_bulk/make_input.py')
mi = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mi)

FIELDS = ('file_id,file_name,file_size,md5sum,access,data_type,analysis.workflow_type,cases.submitter_id,'
          'cases.project.project_id,cases.disease_type,cases.primary_site,cases.samples.sample_type,'
          'cases.samples.submitter_id,cases.samples.tissue_type')


def f_in(field, vals):
    return {'op': 'in', 'content': {'field': field, 'value': vals}}


def files(filters):
    body = {'filters': filters, 'fields': FIELDS, 'size': 20000, 'format': 'JSON'}
    r = urllib.request.Request('https://api.gdc.cancer.gov/files', json.dumps(body).encode(),
                               {'Content-Type': 'application/json'})
    d = json.load(urllib.request.urlopen(r))['data']
    assert d['pagination']['total'] <= 20000, d['pagination']
    return d['hits']


def rows(hits, kind):
    out = []
    for h in hits:
        assert len(h['cases']) == 1, h['file_id']
        c = h['cases'][0]
        s = c.get('samples', [])
        out.append(dict(kind=kind, file_id=h['file_id'], file_name=h['file_name'], md5sum=h['md5sum'],
                        file_size=h['file_size'], access=h['access'], case=c['submitter_id'],
                        primary_site=c.get('primary_site', ''), disease_type=c.get('disease_type', ''),
                        sample_ids=';'.join(x['submitter_id'] for x in s),
                        sample_types=';'.join(x.get('sample_type', '') for x in s),
                        tumour_samples=';'.join(x['submitter_id'] for x in s if x.get('tissue_type') == 'Tumor')))
    return out


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    proj = f_in('cases.project.project_id', ['CPTAC-3'])
    rna = files({'op': 'and', 'content': [proj, f_in('data_type', ['Gene Expression Quantification']),
                                          f_in('analysis.workflow_type', ['STAR - Counts']),
                                          f_in('cases.samples.sample_type', ['Primary Tumor', 'Solid Tissue Normal'])]})
    asc = files({'op': 'and', 'content': [proj, f_in('data_type', ['Allele-specific Copy Number Segment']),
                                          f_in('analysis.workflow_type', ['AscatNGS'])]})
    t = pd.DataFrame(rows(rna, 'rna') + rows(asc, 'ascat'))
    t['cohort'] = [mi.cohort(s, d) for s, d in zip(t.primary_site, t.disease_type)]
    t = t[t.cohort.notna()].sort_values(['kind', 'file_id'])
    t.to_csv(out / 'files.tsv', sep='\t', index=False)
    st = json.load(urllib.request.urlopen('https://api.gdc.cancer.gov/status'))
    (out / 'gdc_status.json').write_text(json.dumps(st, indent=1))
    print(t.groupby(['cohort', 'kind']).size().unstack(), '\nGB', t.file_size.sum() / 1e9, '\naccess', t.access.unique())


if __name__ == '__main__':
    main(sys.argv[1])
