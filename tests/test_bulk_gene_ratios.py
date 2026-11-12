"""Direct gene exports retain focal observations and literal axis identities."""
import json

import numpy as np
import pandas as pd
import pytest

from haplodose.bulk import _direct_gene_ratio, run_bulk


def test_direct_gene_effect_retains_contrast_before_nested_medians():
    # Two triplets form bins. Their medians erase the nonzero focal genes.
    counts=np.array([0.,0.,7.,0.,0.,15.])
    expected=np.zeros(6)
    result,centre=_direct_gene_ratio(counts,expected,[np.arange(3),np.arange(3,6)])
    independent=np.array([0.,0.,np.log(15.),0.,0.,np.log(31.)])
    np.testing.assert_array_equal(result,independent)
    assert centre==0.
    assert np.median(result[:3])==np.median(result[3:])==0.
    assert np.median(result[[2,5]])>0.


@pytest.mark.parametrize('bin_stat',['sum','median'])
def test_bulk_gene_export_reference_and_identity(tmp_path,bin_stat):
    genes=['G3','G0','G5','G1','G4','G2']
    sample_ids=['T','001','N2','N1']
    values=np.array([[10,4,5,7],[20,8,3,5],[30,0,0,4],
                     [12,5,9,3],[18,6,4,8],[2,4,7,9]])
    counts=pd.DataFrame(values,index=genes,columns=sample_ids)
    counts.to_csv(tmp_path/'counts.tsv',sep='\t')
    pd.DataFrame(dict(sample=sample_ids,role=['tumour','normal','normal','normal'])).to_csv(tmp_path/'samples.tsv',sep='\t',index=False)
    pd.DataFrame(dict(gene=['G'+str(i) for i in range(6)],gene_symbol=['G'+str(i) for i in range(6)],chromosome='chr1',
                      start=np.arange(6)*1000+130000000,end=np.arange(6)*1000+130000500)).to_csv(tmp_path/'positions.tsv',sep='\t',index=False)
    with pytest.warns(UserWarning,match='No panel size establishes calibration'):
        out=run_bulk(tmp_path/'counts.tsv',tmp_path/'samples.tsv',tmp_path/'positions.tsv',tmp_path/'out',
                     bin_genes=2,min_cpm=0,n_programmes=0,bin_stat=bin_stat,counts_kind='raw',
                     genome_build='hg38',coordinate_system='zero-based-half-open',gene_namespace='custom')
    mapping=pd.read_csv(out/'genes.tsv.gz',sep='\t')
    observed=np.load(out/'gene_rna_log_ratio.npy')
    exported=pd.read_csv(out/'samples.tsv',sep='\t',dtype={'sample':str})['sample'].tolist()
    assert exported==sample_ids
    metadata=pd.read_csv(out/'gene_ratio_metadata.tsv',sep='\t',dtype={'sample':str})
    assert metadata['sample'].tolist()==sample_ids
    for row,sample in enumerate(sample_ids):
        references=[name for name in sample_ids[1:] if name!=sample]
        total=sum(int(counts[name].sum()) for name in references)
        library=int(counts[sample].sum())
        uncentred=[]
        for gene in mapping.gene:
            ref=sum(int(counts.loc[gene,name]) for name in references)/total
            uncentred.append(np.log((int(counts.loc[gene,sample])+.5)/(library*ref+.5)))
        uncentred=np.array(uncentred)
        centre=np.median([np.median(uncentred[mapping.bin==b]) for b in sorted(mapping.bin.unique())])
        np.testing.assert_array_equal(observed[row],uncentred-centre)
        assert metadata.loc[row,'n_reference_samples']==len(references)
        assert metadata.loc[row,'reference_input_count']==total
    run=json.loads((out/'run.json').read_text())['gene_rna_log_ratio']
    assert run['calibration_status']=='unvalidated'
    assert run['programme_correction'] is False and run['smoothing'] is False
    assert run['feature_selection']=='common_all_normal_feature_universe_not_leave_one_out'
    manifest=json.loads((out/'run_manifest.json').read_text())
    assert 'gene_rna_log_ratio.npy' in manifest['output']['files']
