"""REBC-THYR bulk sub-arm benchmark (prespecified claim 7, docs/STATUS.md): benchmarks/cptac_bulk/score.py unchanged
(same scorer, primary test, Holm over 2, bootstrap seed). Afterwards the labels in primary.json are corrected to
REBC-THYR and the claim 7 power floor (>= 20 tumours with a focal-loss window) is applied.
Usage: score.py OUT INPUT_ROOT ICNV_ROOT
"""
import importlib.util
import json
import sys
from pathlib import Path

H = Path('/gpfs/commons/home/jameslee/Haplodose')
_spec = importlib.util.spec_from_file_location('cptac_score', H / 'benchmarks/cptac_bulk/score.py')
cs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cs)
MIN_TUMOURS = 20


def main(out, inp, icnv):
    cs.main(out, inp, icnv)
    p = Path(out) / 'primary.json'
    res = json.loads(p.read_text())
    res['unit'] = 'tumour (REBC-THYR case), tumours with >= 1 focal-loss window'
    res['truth_note'] = 'GDC REBC-THYR AscatNGS segments (GRCh38); subarm GTF37 string inherited from the CPTAC driver'
    a = res['all']
    n = min(c['n_lines'] for c in a['comparisons'])   # tumours with >= 1 focal-loss window and paired recalls
    res['power_floor'] = dict(min_tumours=MIN_TUMOURS, n=n, evaluable=n is not None and n >= MIN_TUMOURS)
    res['claim7_holds'] = bool(res['power_floor']['evaluable'] and a.get('claim_holds'))
    p.write_text(json.dumps(res, indent=1))
    print('CLAIM7', json.dumps(res['power_floor']), 'holds', res['claim7_holds'])


if __name__ == '__main__':
    main(*sys.argv[1:4])
