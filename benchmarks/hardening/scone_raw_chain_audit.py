"""Independent reconstruction of one staged library's bounded raw-prefix report."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

from haplodose.provenance import start_run,finish_run,file_record

MOTIFS={'TCATG':'DNA','AACTT':'RNA_author_FASTA','ATCGT':'RNA_GEO_protocol',
        'CATGA':'DNA_reverse_complement','AAGTT':'RNA_FASTA_reverse_complement',
        'ACGAT':'RNA_GEO_reverse_complement'}


def increment(table,key):
    table[key]=table.get(key,0)+1


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True);p.add_argument('--metadata-task',type=Path,required=True)
    p.add_argument('--library-class',choices=['legacy_whole','later_nucleus'],default='legacy_whole')
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    directory=a.task/('chain_v1' if a.library_class=='legacy_whole' else 'chain_nucleus_v1')
    design='AUDIT_DESIGN.md' if a.library_class=='legacy_whole' else 'NUCLEUS_DESIGN.md'
    run=start_run(a.out,dict(script=Path(__file__),design=a.task/design,
                            source_run=directory/'run_manifest.json',
                            original_ENA=a.metadata_task/'metadata_v1/ENA_runs.tsv'),dict(library_class=a.library_class),
                  dict(scope='Independent raw-prefix integer counts; no alignments or CNA inference'),mode='scone_raw_chain_audit')
    try:
        manifest=json.loads((directory/'run_manifest.json').read_text());assert manifest['status']=='complete'
        for kind in ('input','output'):
            for record in manifest[kind]['files'].values():
                assert file_record(record['path'])['sha256']==record['sha256']
        inventory=json.loads((directory/'output_files.json').read_text())
        for record in inventory.values():assert file_record(record['path'])['sha256']==record['sha256']
        with (a.metadata_task/'metadata_v1/ENA_runs.tsv').open() as handle:
            rows=list(csv.DictReader(handle,delimiter='\t'))
        pattern=r': HCT116 scONEseq \d+ Homo sapiens' if a.library_class=='legacy_whole' else r': HCT116, extracted nuclei, cell \d+ Homo sapiens'
        eligible=[r for r in rows if re.search(pattern,r['experiment_title'])
                  and r['library_strategy']=='OTHER']
        selected=sorted(eligible,key=lambda r:(r['sample_alias'],r['run_accession']))[0]
        report=json.loads((directory/'report.json').read_text())
        expected=('GSM5609521','SRR16194948') if a.library_class=='legacy_whole' else ('GSM6683025','SRR22047451')
        assert report['selected_sample']==selected['sample_alias']==expected[0]
        assert report['selected_run']==selected['run_accession']==expected[1]
        assert 0<report['pairs_examined']<=250000==report['maximum_pairs']
        sizes=list(map(int,selected['fastq_bytes'].split(';')));hashes=selected['fastq_md5'].split(';')
        for record,size,expected in zip(report['files'],sizes,hashes):
            path=Path(record['path']);assert path.stat().st_size==size==record['expected_bytes']
            digest=hashlib.md5()
            with path.open('rb') as handle:
                for block in iter(lambda:handle.read(4*1024*1024),b''):digest.update(block)
            assert digest.hexdigest()==expected==record['md5']
        prefix=[{},{}];lengths=[{},{}];offsets=[{name:{} for name in MOTIFS.values()} for _ in range(2)]
        pairs={};n=report['pairs_examined']
        with gzip.open(report['files'][0]['path'],'rt',encoding='ascii') as first,gzip.open(report['files'][1]['path'],'rt',encoding='ascii') as second:
            for i in range(n):
                names=[];categories=[]
                for mate,handle in enumerate((first,second)):
                    header=handle.readline().rstrip('\r\n');sequence=handle.readline().rstrip('\r\n').upper()
                    plus=handle.readline();quality=handle.readline().rstrip('\r\n')
                    assert header.startswith('@') and plus.startswith('+') and sequence and len(sequence)==len(quality)
                    names.append(re.sub(r'/[12]$','',header[1:].split()[0]))
                    increment(prefix[mate],sequence[:5]);increment(lengths[mate],str(len(sequence)))
                    categories.append(MOTIFS.get(sequence[:5],'other_prefix'))
                    # Independent per-position slices, instead of str.find.
                    for position in range(min(31,max(0,len(sequence)-4))):
                        motif=sequence[position:position+5]
                        if motif in MOTIFS:increment(offsets[mate][MOTIFS[motif]],str(position))
                assert names[0]==names[1]
                if i<3:assert names[0]==report['first_three_pair_ids'][i]
                increment(pairs,'|'.join(categories))
        assert prefix==report['prefix_5bp_counts']
        assert lengths==report['read_length_counts']
        assert offsets==report['motif_occurrences_start_0_through_30']
        assert pairs==report['paired_prefix_category_counts']
        metadata=(directory/'sample.soft.txt').read_text()
        titles=re.findall(r'^!Sample_title = (.+)$',metadata,re.M)
        descriptions=re.findall(r'^!Sample_description = (.+)$',metadata,re.M)
        (a.out/'report.json').write_text(json.dumps(dict(selected_library_exact=True,
                    archive_MD5_and_sizes_verified=True,recursive_hashes_verified=True,pairs_reconstructed=n,
                    all_prefix_length_offset_and_pair_counts_exact=True,titles=titles,descriptions=descriptions,
                    literal_GSM_and_SRR_are_identity_keys=True,numeric_title_description_mapping='unverified',
                    biological_null_specimens_verified=0,protocol_chosen=False,alignment_or_CNA_inference=False),indent=2)+'\n')
        finish_run(a.out,run,backend='cpu')
    except Exception as error:
        finish_run(a.out,run,error=error);raise


if __name__=='__main__':main()
