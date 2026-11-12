"""Prepare the two fixed bounded library subsets with explicit candidate layouts."""
import argparse
from collections import Counter
from contextlib import ExitStack
import gzip
import json
from pathlib import Path

from haplodose.provenance import start_run, finish_run, file_record

ME=b'AGATGTGTATAAGAGACAG'
LAYOUTS=[dict(sample='GSM5609521',run='SRR16194948',chain='chain_v1',rna_tag='ATCGT',umi=5),
         dict(sample='GSM6683025',run='SRR22047451',chain='chain_nucleus_v1',rna_tag='AACTT',umi=8)]


def read(handle):
    header=handle.readline();sequence=handle.readline().rstrip(b'\r\n')
    plus=handle.readline();quality=handle.readline().rstrip(b'\r\n')
    assert header.startswith(b'@') and plus.startswith(b'+') and len(sequence)==len(quality)>0
    name=header[1:].split()[0]
    if name.endswith((b'/1',b'/2')):name=name[:-2]
    return name,sequence,quality


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--raw-task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    inputs=dict(script=Path(__file__),design=args.task/'PREPARE_DESIGN.md')
    for name in ('chain_v1','chain_nucleus_v1','audit_v1','audit_nucleus_v1','layout_v1','layout_nucleus_v1','rna_layout_v1'):
        inputs[name]=args.raw_task/name/'run_manifest.json'
    run=start_run(args.out,inputs,dict(schema_version=1,layouts=LAYOUTS,pairs_per_library=250000),
                  dict(scope='Candidate read layers only, not molecule counts, genomic truth or calibration'),mode='scone_prepare_reads')
    try:
        for name,path in inputs.items():
            if name in ('script','design'):continue
            manifest=json.loads(path.read_text());assert manifest['status']=='complete'
            for record in manifest['output']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        reports={}
        for layout in LAYOUTS:
            source=json.loads((args.raw_task/layout['chain']/'report.json').read_text())
            assert source['selected_sample']==layout['sample'] and source['selected_run']==layout['run']
            assert source['pairs_examined']==250000
            folder=args.out/layout['run'];folder.mkdir();counts=Counter();seen=set();prefix=Counter()
            with ExitStack() as stack:
                raw=[stack.enter_context(gzip.open(record['path'],'rb')) for record in source['files']]
                outputs={}
                for layer in ('RNA','DNA'):
                    for mate in (1,2):
                        handle=stack.enter_context((folder/f'{layer}_R{mate}.fastq.gz').open('xb'))
                        outputs[layer,mate]=stack.enter_context(gzip.GzipFile(filename='',mode='wb',fileobj=handle,mtime=0))
                for _ in range(250000):
                    first,second=[read(handle) for handle in raw]
                    assert first[0]==second[0] and first[0] not in seen;seen.add(first[0])
                    name,sequence,quality=second;tag=sequence[:5].decode('ascii');prefix[tag]+=1
                    if tag==layout['rna_tag']:layer='RNA'
                    elif tag=='TCATG':layer='DNA'
                    else:counts['other_prefix']+=1;continue
                    stop=5+layout['umi'];umi=sequence[5:stop]
                    if len(umi)!=layout['umi'] or any(base not in b'ACGT' for base in umi):
                        counts['invalid_candidate_UMI']+=1;continue
                    if layer=='DNA':
                        if sequence[stop:stop+len(ME)]!=ME:counts['DNA_missing_fixed_ME']+=1;continue
                        stop+=len(ME)
                    if not sequence[stop:]:counts['empty_post_trim']+=1;continue
                    emitted_name=name+b':umi_'+umi
                    for mate,seq,qual in [(1,first[1],first[2]),(2,sequence[stop:],quality[stop:])]:
                        outputs[layer,mate].write(b'@'+emitted_name+b'\n'+seq+b'\n+\n'+qual+b'\n')
                    counts[layer]+=1
            assert sum(counts.values())==len(seen)==250000
            assert dict(prefix)==source['prefix_5bp_counts'][1]
            reports[layout['run']]=dict(layout=layout,counts=dict(counts),unique_input_pairs=len(seen),
                                       files={f.name:file_record(f) for f in sorted(folder.iterdir())})
        (args.out/'report.json').write_text(json.dumps(dict(schema_version=1,libraries=reports,
                 molecule_counts=False,alignment_run=False,biological_null_specimens_verified=0),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
