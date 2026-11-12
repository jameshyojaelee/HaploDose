"""Verify the unchanged hard-control workflow against the prior frozen thin replay."""
import argparse,json
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
p.add_argument('--current',type=Path,required=True);p.add_argument('--previous',type=Path,required=True)
a=p.parse_args();a.out.mkdir(exist_ok=False);count=0;groups=[]
for role,n in [('cal_nb',199),('cal_state',199),('eval_nb',200),('eval_state',200),('eval_two',200)]:
 old=[json.loads(x) for x in (a.previous/role/'records.jsonl').read_text().splitlines()]
 new=[json.loads(x) for x in (a.current/'hard'/role/'records.jsonl').read_text().splitlines()]
 if len(old)!=n or len(new)!=n:raise ValueError('Incomplete control replay')
 for before,after in zip(old,new):
    for name in ('condition','role','replicate','seed_key','input_sha256'):
        if before[name]!=after[name]:raise ValueError('Control input identity differs: '+str(after['seed_key']))
    ob=np.load(before['assignments']['path']);nb=np.load(after['assignments']['path'])
    for name in ('simulated_CNA_labels','raw_MAP','raw_R','latent_state','eligible_tumour_groups'):
        if not np.array_equal(ob[name],nb[name]):raise ValueError('Control array changed: '+name+' '+str(after['seed_key']))
    if before['original_final_difference_score']!=after['score']:
        raise ValueError('Control score changed: '+str(after['seed_key']))
    count+=1
 groups.append(dict(role=role,identical_specimens=n))
report=dict(identical_input_probability_assignment_score_specimens=count,roles=groups,
            scope='Exact preservation of prior thin-mode raw workflow; prior event-observer score is not reused as final-group score')
(a.out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
