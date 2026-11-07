"""Version-2 retrospective REBC scoring. The read-only original scorer is preserved."""
import importlib.util
import json
import sys
from pathlib import Path

H = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('cptac_score_v2', H / 'benchmarks/cptac_bulk/score.py')
cs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cs)
MIN_TUMOURS = 20


def main(out, inp, icnv):
    cs.main(out, inp, icnv)
    path = Path(out) / 'primary.json'
    res = json.loads(path.read_text())
    res['unit'] = 'tumour (REBC-THYR case), tumours with >= 1 callable focal-loss window'
    res['truth_note'] = 'GDC REBC-THYR AscatNGS segments (GRCh38)'
    n = min(c['n_lines'] for c in res['all']['comparisons'])
    res['power_floor'] = dict(min_tumours=MIN_TUMOURS, n=n, evaluable=n >= MIN_TUMOURS)
    res['descriptive_rule_met'] = bool(n >= MIN_TUMOURS and res['all']['descriptive_rule_met'])
    res['claim7_holds'] = False
    path.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main(*sys.argv[1:4])
