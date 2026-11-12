"""Bounded source-specified RNA-primer signatures, with independent exact matching."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import re

from haplodose.provenance import start_run, finish_run, file_record

RULES={'GGGHN':('GGG',set('ACT')),'TTTVN':('TTT',set('ACG'))}
REGEX={name:re.compile('^'+triplet+'['+''.join(sorted(bases))+'][ACGT]')
       for name,(triplet,bases) in RULES.items()}


def signatures(sequence,offset):
    tail=sequence[offset:]
    observed={name:len(tail)>=5 and tail[:3]==triplet and tail[3] in bases and tail[4] in 'ACGT'
              for name,(triplet,bases) in RULES.items()}
    n_t=0
    for base in tail:
        if base!='T':break
        n_t+=1
    observed['polyT_VN']=n_t>=20 and len(tail)>=n_t+2 and tail[n_t] in 'ACG' and tail[n_t+1] in 'ACGT'
    return observed


def independent(sequence,offset):
    tail=sequence[offset:]
    return {**{name:bool(pattern.match(tail)) for name,pattern in REGEX.items()},
            'polyT_VN':bool(re.match(r'^T{20,}[ACG][ACGT]',tail))}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    inputs=dict(script=Path(__file__),design=args.task/'RNA_LAYOUT_DESIGN.md',
                primary_table=args.task/'primer_cloud_v1/run_manifest.json')
    for name in ('chain_v1','audit_v1','chain_nucleus_v1','audit_nucleus_v1'):
        inputs[name]=args.task/name/'run_manifest.json'
    run=start_run(args.out,inputs,dict(candidates=[5,8],maximum_pairs=250000,triplet_offsets=list(range(41))),
                  dict(scope='RNA primer-layout development diagnostic only; no molecule or CNA inference'),mode='scone_rna_layout')
    try:
        for name,path in inputs.items():
            if name in ('script','design'):continue
            manifest=json.loads(path.read_text());assert manifest['status']=='complete'
            for record in manifest['output']['files'].values():
                assert file_record(record['path'])['sha256']==record['sha256']
        controls=0
        for umi in ('ACGTA','ACGTACGT'):
            for flank,target in [('GGGAC','GGGHN'),('TTTAC','TTTVN'),('T'*20+'AC','polyT_VN')]:
                sequence='AACTT'+umi+flank+'CGACAC'
                expected={name:name==target for name in (*RULES,'polyT_VN')}
                assert signatures(sequence,5+len(umi))==independent(sequence,5+len(umi))==expected
                controls+=1
        for flank in ('GGGGA','TTTTA','GGGAN'):
            assert not any(signatures(flank,0).values())
            assert signatures(flank,0)==independent(flank,0);controls+=1
        results={}
        for label,chain in [('legacy','chain_v1'),('nucleus','chain_nucleus_v1')]:
            source=json.loads((args.task/chain/'report.json').read_text());n=source['pairs_examined'];assert n==250000
            record=source['files'][1];assert file_record(record['path'])['sha256']==record['sha256']
            strata={tag:dict(reads=0,short_for_candidate=Counter(),candidate_counts={str(umi):Counter() for umi in (5,8)},
                           triplet_offsets={triplet:Counter() for triplet in ('GGG','TTT')}) for tag in ('ATCGT','AACTT','TCATG')}
            with gzip.open(record['path'],'rt',encoding='ascii') as handle:
                for _ in range(n):
                    header=handle.readline();sequence=handle.readline().rstrip('\r\n').upper()
                    plus=handle.readline();quality=handle.readline().rstrip('\r\n')
                    assert header.startswith('@') and plus.startswith('+') and len(sequence)==len(quality)>0
                    tag=sequence[:5]
                    if tag not in strata:continue
                    row=strata[tag];row['reads']+=1
                    for umi in (5,8):
                        observed=signatures(sequence,5+umi)
                        assert observed==independent(sequence,5+umi)
                        row['short_for_candidate'][str(umi)]+=len(sequence)<5+umi+5
                        for name,present in observed.items():row['candidate_counts'][str(umi)][name]+=int(present)
                    for triplet in ('GGG','TTT'):
                        direct=Counter(i for i in range(41) if sequence[i:i+3]==triplet)
                        regex=Counter(m.start() for m in re.finditer('(?='+triplet+')',sequence) if m.start()<=40)
                        assert direct==regex
                        row['triplet_offsets'][triplet].update({str(k):v for k,v in direct.items()})
            for tag,row in strata.items():assert row['reads']==source['prefix_5bp_counts'][1].get(tag,0)
            results[label]=dict(pairs_examined=n,strata=strata)
        (args.out/'report.json').write_text(json.dumps(dict(libraries=results,synthetic_controls=controls,
                    independent_actual_read_patterns_exact=True,biological_null_specimens_verified=0,
                    processing_protocol_chosen=False,alignment_or_quantification=False),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
