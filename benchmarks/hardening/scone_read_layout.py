"""Bounded, source-specified adapter-offset check; no molecule counting or alignment."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import re

from haplodose.provenance import start_run,finish_run,file_record

ME='AGATGTGTATAAGAGACAG'
REVERSE_ME='CTGTCTCTTATACACATCT'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True)
    p.add_argument('--library-class',choices=['legacy_whole','later_nucleus'],default='legacy_whole')
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    chain='chain_v1' if a.library_class=='legacy_whole' else 'chain_nucleus_v1'
    audit='audit_v1' if a.library_class=='legacy_whole' else 'audit_nucleus_v1'
    design='LAYOUT_DESIGN.md' if a.library_class=='legacy_whole' else 'NUCLEUS_DESIGN.md'
    run=start_run(a.out,dict(script=Path(__file__),design=a.task/design,
                            source_run=a.task/chain/'run_manifest.json',
                            audit=a.task/audit/'run_manifest.json',
                            archived_source=a.task/'archived_source_v1/run_manifest.json'),
                  dict(library_class=a.library_class,candidate_UMI_lengths=[5,8],tag_length=5,ME=ME,reverse_ME=REVERSE_ME,limit_pairs=250000),
                  dict(scope='One development-library layout diagnostic, not quantified RNA/DNA truth'),mode='scone_read_layout')
    try:
        for directory in (chain,audit,'archived_source_v1'):
            m=json.loads((a.task/directory/'run_manifest.json').read_text());assert m['status']=='complete'
            for record in m['output']['files'].values():
                assert file_record(record['path'])['sha256']==record['sha256']
        for umi in ('ACGTA','ACGTACGT'):
            seq='TCATG'+umi+ME+'GATTACA'
            observed=[i for i in range(41) if seq[i:i+len(ME)]==ME]
            assert observed==[5+len(umi)]
            assert observed==[m.start() for m in re.finditer('(?='+ME+')',seq)]
        source=json.loads((a.task/chain/'report.json').read_text());n=source['pairs_examined']
        path=Path(source['files'][1]['path']);assert file_record(path)['sha256']==source['files'][1]['sha256']
        names=(('DNA_tag','RNA_GEO_tag','other_prefix') if a.library_class=='legacy_whole'
               else ('DNA_tag','RNA_GEO_tag','RNA_author_tag','other_prefix'))
        positions={name:[Counter() for _ in range(60)] for name in names}
        offsets={name:{'ME':Counter(),'reverse_ME':Counter()} for name in names}
        sizes=Counter();independent={name:{'ME':Counter(),'reverse_ME':Counter()} for name in names}
        with gzip.open(path,'rt',encoding='ascii') as handle:
            for _ in range(n):
                header=handle.readline();sequence=handle.readline().rstrip('\r\n').upper()
                plus=handle.readline();quality=handle.readline().rstrip('\r\n')
                assert header.startswith('@') and plus.startswith('+') and len(sequence)==len(quality)>0
                name='DNA_tag' if sequence.startswith('TCATG') else 'RNA_GEO_tag' if sequence.startswith('ATCGT') else 'other_prefix'
                if a.library_class=='later_nucleus' and sequence.startswith('AACTT'):
                    name='RNA_author_tag'
                sizes[name]+=1
                for i,base in enumerate(sequence[:60]):positions[name][i][base]+=1
                for kind,motif in (('ME',ME),('reverse_ME',REVERSE_ME)):
                    for i in range(41):
                        if sequence[i:i+len(motif)]==motif:offsets[name][kind][i]+=1
                    # Different matching mechanism for the actual reads, not
                    # just the two synthetic offset fixtures.
                    for match in re.finditer('(?='+motif+')',sequence):
                        if match.start()<=40:independent[name][kind][match.start()]+=1
        assert offsets==independent
        assert sizes['DNA_tag']==source['prefix_5bp_counts'][1]['TCATG']
        assert sizes['RNA_GEO_tag']==source['prefix_5bp_counts'][1]['ATCGT']
        if a.library_class=='later_nucleus':
            assert sizes['RNA_author_tag']==source['prefix_5bp_counts'][1].get('AACTT',0)
        assert sum(sizes.values())==n
        result=dict(pairs_examined=n,read_mate=2,tag_stratum_counts=dict(sizes),
                    positional_base_counts={name:[dict(sorted(c.items())) for c in rows] for name,rows in positions.items()},
                    exact_adapter_offsets={name:{kind:dict(sorted(c.items())) for kind,c in kinds.items()} for name,kinds in offsets.items()},
                    candidate_DNA_layout_counts={str(umi):offsets['DNA_tag']['ME'][5+umi] for umi in (5,8)},
                    independent_regex_counts_exact=True,synthetic_offset_controls=2,
                    interpretation='Fixed adapter offset can distinguish candidate DNA layouts; it does not prove RNA UMI identity or an executed workflow',
                    biological_null_specimens_verified=0,protocol_chosen=False,alignment_or_quantification=False)
        (a.out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
        finish_run(a.out,run,backend='cpu')
    except Exception as error:
        finish_run(a.out,run,error=error);raise


if __name__=='__main__':main()
