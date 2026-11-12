"""Independent raw-read regex reconstruction of every prepared paired record."""
import argparse
from collections import Counter
from contextlib import ExitStack
import gzip
import json
from pathlib import Path
import re

from haplodose.provenance import start_run, finish_run, file_record


def record(handle):
    header=handle.readline().rstrip('\r\n');sequence=handle.readline().rstrip('\r\n')
    plus=handle.readline().rstrip('\r\n');quality=handle.readline().rstrip('\r\n')
    assert header.startswith('@') and plus.startswith('+') and len(sequence)==len(quality)>0
    return re.sub(r'/[12]$','',header[1:].split()[0]),sequence,quality


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--raw-task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False);prepared=args.task/'prepare_v1'
    run=start_run(args.out,dict(script=Path(__file__),design=args.task/'PREPARE_DESIGN.md',
                            prepared_run=prepared/'run_manifest.json'),{},
                  dict(scope='Exact conservation of prepared sequences and qualities, not assay validation'),mode='scone_prepare_reads_audit')
    try:
        manifest=json.loads((prepared/'run_manifest.json').read_text());assert manifest['status']=='complete'
        for section in ('input','output'):
            for item in manifest[section]['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
        report=json.loads((prepared/'report.json').read_text());verified={}
        # Literal independent constants, not the producer's chosen layout fields.
        for accession,sample,chain,tag,n_umi in [('SRR16194948','GSM5609521','chain_v1','ATCGT',5),
                                               ('SRR22047451','GSM6683025','chain_nucleus_v1','AACTT',8)]:
            emitted=report['libraries'][accession]
            assert emitted['layout']==dict(sample=sample,run=accession,chain=chain,rna_tag=tag,umi=n_umi)
            source=json.loads((args.raw_task/chain/'report.json').read_text());assert source['selected_run']==accession
            for item in source['files']:
                assert file_record(item['path'])['sha256']==item['sha256']
            for item in emitted['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
            patterns={layer:re.compile('^'+prefix+'(?P<umi>[ACGT]{'+str(n_umi)+'})'+flank+'(?P<insert>.+)$')
                      for layer,prefix,flank in [('RNA',tag,''),('DNA','TCATG','AGATGTGTATAAGAGACAG')]}
            counts=Counter();ids=set()
            with ExitStack() as stack:
                raw=[stack.enter_context(gzip.open(item['path'],'rt',encoding='ascii')) for item in source['files']]
                output={(layer,mate):stack.enter_context(gzip.open(prepared/accession/f'{layer}_R{mate}.fastq.gz','rt',encoding='ascii'))
                        for layer in ('RNA','DNA') for mate in (1,2)}
                for _ in range(250000):
                    first,second=[record(handle) for handle in raw]
                    assert first[0]==second[0] and first[0] not in ids;ids.add(first[0])
                    sequence=second[1];matches=[(layer,pattern.match(sequence)) for layer,pattern in patterns.items()]
                    matches=[(layer,match) for layer,match in matches if match is not None]
                    assert len(matches)<=1
                    if matches:
                        layer,match=matches[0];name=first[0]+':umi_'+match.group('umi')
                        assert record(output[layer,1])==(name,first[1],first[2])
                        assert record(output[layer,2])==(name,match.group('insert'),second[2][match.start('insert'):])
                        counts[layer]+=1
                    elif sequence[:5] not in (tag,'TCATG'):counts['other_prefix']+=1
                    elif not re.match(r'^[ACGT]{'+str(n_umi)+'}$',sequence[5:5+n_umi]):counts['invalid_candidate_UMI']+=1
                    elif sequence[:5]=='TCATG' and not re.match('^TCATG[ACGT]{'+str(n_umi)+'}AGATGTGTATAAGAGACAG',sequence):
                        counts['DNA_missing_fixed_ME']+=1
                    else:counts['empty_post_trim']+=1
                assert all(handle.readline()=='' for handle in output.values()),'Unexpected leftover prepared record'
            assert dict(counts)==emitted['counts'] and sum(counts.values())==len(ids)==emitted['unique_input_pairs']==250000
            verified[accession]=dict(counts=dict(counts),all_emitted_records_exact=True,unique_pairs=len(ids))
        (args.out/'report.json').write_text(json.dumps(dict(libraries=verified,all_raw_and_output_hashes_verified=True,
                 biological_null_specimens_verified=0,molecule_counts=False),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
