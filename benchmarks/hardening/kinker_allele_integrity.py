"""Read-only donor-restricted pileup and exported-phase integrity trace."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from haplodose.provenance import start_run, finish_run


def vcf(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as stream:
        for line in stream:
            if line.startswith('#CHROM'):
                names = line.strip().lstrip('#').split('\t')
                return pd.read_csv(stream, sep='\t', names=names, dtype=str)
    raise ValueError(f'VCF header absent: {path}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args(); a.out.mkdir(exist_ok=False)
    selected = json.loads((a.task/'selection.json').read_text())['selected']
    directory = Path(selected['path']); identifier = selected['identifier']
    root = a.task.parents[1]
    inputs = dict(script=Path(__file__), selection=a.task/'selection.json', design=a.task/'DESIGN.md',
                  phase_code=root/'benchmarks/kinker/alleles/phase_pool.R',
                  exported=directory/(identifier+'_allele_counts.tsv.gz'))
    channels = sorted((directory/'pileup').iterdir())
    for channel in channels:
        original = root/'results/kinker_numbat_20261102/pileup'/channel.name
        for name in ['cell_map.tsv', 'cellSNP.samples.tsv', 'cellSNP.base.vcf',
                     'cellSNP.tag.AD.mtx', 'cellSNP.tag.DP.mtx', 'cellSNP.tag.OTH.mtx']:
            inputs[channel.name+':original:'+name] = original/name
        for name in ['cellSNP.samples.tsv', 'cellSNP.base.vcf', 'cellSNP.tag.AD.mtx', 'cellSNP.tag.DP.mtx']:
            inputs[channel.name+':subset:'+name] = channel/name
    for chrom in range(1,23):
        for suffix in ['.vcf.gz', '.phased.vcf.gz']:
            inputs[f'{chrom}{suffix}'] = directory/'phasing'/f'{identifier}_chr{chrom}{suffix}'
    run = start_run(a.out, inputs, dict(identifier=identifier, channels=[p.name for p in channels]),
                    dict(scope='One historical donor integrity trace; no biological validation'),
                    mode='kinker_allele_integrity')
    try:
        exported = pd.read_csv(inputs['exported'], sep='\t', dtype={'CHROM':str,'GT':str})
        checks = []; pileup = {}
        for channel in channels:
            original = root/'results/kinker_numbat_20261102/pileup'/channel.name
            mapping = pd.read_csv(original/'cell_map.tsv', sep='\t')
            mapping = mapping[mapping.line == identifier]
            barcodes = pd.read_csv(original/'cellSNP.samples.tsv', header=None)[0]
            assert barcodes.is_unique and mapping.bam_barcode.is_unique and mapping.cell_id.is_unique
            indices = pd.Index(barcodes).get_indexer(mapping.bam_barcode)
            missing = mapping.loc[indices < 0, 'bam_barcode'].tolist()
            assert not missing, dict(channel=channel.name, missing=missing)
            matrices = {k:mmread(original/f'cellSNP.tag.{k}.mtx').tocsr()[:,indices] for k in ['AD','DP','OTH']}
            totals = {k:np.asarray(v.sum(axis=1)).ravel() for k,v in matrices.items()}
            keep = totals['DP'] + totals['OTH'] >= 2
            full = vcf(original/'cellSNP.base.vcf'); sub = vcf(channel/'cellSNP.base.vcf')
            assert len(full) == len(keep)
            pd.testing.assert_frame_equal(full.loc[keep,['CHROM','POS','REF','ALT']].reset_index(drop=True),
                                          sub[['CHROM','POS','REF','ALT']])
            cells = pd.read_csv(channel/'cellSNP.samples.tsv', header=None)[0]
            assert cells.tolist() == mapping.cell_id.tolist()
            for k in ['AD','DP']:
                saved = mmread(channel/f'cellSNP.tag.{k}.mtx').tocsr()
                assert saved.shape == matrices[k][keep].shape and (saved != matrices[k][keep]).nnz == 0
            expected_info = ['AD='+str(int(ad))+';DP='+str(int(dp))+';OTH='+str(int(oth))
                             for ad,dp,oth in zip(totals['AD'][keep], totals['DP'][keep], totals['OTH'][keep])]
            assert sub.INFO.tolist() == expected_info
            selected_counts = {k:matrices[k][keep].toarray() for k in ['AD','DP']}
            for i, row in sub.iterrows():
                for j, cell in enumerate(cells):
                    pileup[(cell,str(row.CHROM).removeprefix('chr'),int(row.POS),row.REF,row.ALT)] = (
                        int(selected_counts['AD'][i,j]), int(selected_counts['DP'][i,j]))
            checks.append(dict(channel=channel.name, cells=len(cells), retained_snps=int(keep.sum()), missing_barcodes=missing))
        phases = {}; genotype_checks = 0
        for chrom in range(1,23):
            target = vcf(inputs[f'{chrom}.vcf.gz']); phased = vcf(inputs[f'{chrom}.phased.vcf.gz'])
            assert list(target.columns[9:]) == list(phased.columns[9:]) == [identifier]
            def genotypes(frame):
                return {(str(r.CHROM).removeprefix('chr'),int(r.POS),r.REF,r.ALT):str(r[identifier]).split(':')[0]
                        for _,r in frame.iterrows()}
            before, after = genotypes(target), genotypes(phased)
            assert before.keys() == after.keys()
            for key, gt in after.items():
                assert sorted(gt.replace('|','/').split('/')) == sorted(before[key].replace('|','/').split('/'))
                phases[key] = gt; genotype_checks += 1
        for row in exported.itertuples(index=False):
            key = (str(row.CHROM).removeprefix('chr'),int(row.POS),row.REF,row.ALT)
            assert row.GT in ('0|1','1|0') and phases[key] == row.GT
            assert pileup[(row.cell,*key)] == (row.AD,row.DP)
        (a.out/'report.json').write_text(json.dumps(dict(identifier=identifier,channels=checks,
            genotype_phase_checks=genotype_checks, exported_rows_checked=len(exported),
            assembly_provenance='Not established by count/phase equality; inspect original reference metadata',
            executed_historical_source='Current phase script recorded; historical execution identity not independently established',
            biological_validation=False),indent=2))
        finish_run(a.out,run,backend='cpu')
    except Exception as error:
        finish_run(a.out,run,error=error); raise


if __name__ == '__main__':
    main()
