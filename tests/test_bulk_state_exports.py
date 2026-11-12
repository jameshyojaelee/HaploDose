"""Bulk exports preserve conditional per-bin evidence before segment averaging."""
import json

import numpy as np
import pandas as pd
import pytest

from haplodose import bulk
from haplodose.profile import LOG_DOSE


@pytest.mark.parametrize('mode', ['sum', 'median'])
def test_bulk_state_exports_keep_per_bin_fit_and_literal_axes(tmp_path, monkeypatch, mode):
    genes = ['G3','G0','G5','G1','G4','G2']; names = ['T','001','N2','N1']
    pd.DataFrame([[10,4,5,7],[20,8,3,5],[30,2,1,4],[12,5,9,3],[18,6,4,8],[2,4,7,9]],
                 index=genes, columns=names).to_csv(tmp_path/'counts.tsv',sep='\t')
    pd.DataFrame(dict(sample=names,role=['tumour','normal','normal','normal'])).to_csv(tmp_path/'samples.tsv',sep='\t',index=False)
    pd.DataFrame(dict(gene=[f'G{i}' for i in range(6)],gene_symbol=[f'G{i}' for i in range(6)],chromosome='chr1',start=np.arange(6)*1000+130000000,
                      end=np.arange(6)*1000+130000500)).to_csv(tmp_path/'positions.tsv',sep='\t',index=False)
    captured = []
    def fit(y, mu, v, W, reset, params):
        probability = np.array([[.1,.6,.2,.1],[.3,.4,.2,.1],[.05,.75,.1,.1]],np.float32)
        effect = probability @ (.6*LOG_DOSE)
        path = np.ones(3,np.int32)
        captured.append((probability.copy(),effect.copy(),path.copy()))
        return dict(E=y.sum()*mu, shift=np.zeros(3), centre=0., beta=.6, tau=.15,
                    post=probability, dose=effect, path=path)
    monkeypatch.setattr(bulk,'fit_sample',fit)
    with pytest.warns(UserWarning,match='No panel size establishes calibration'):
        out = bulk.run_bulk(tmp_path/'counts.tsv',tmp_path/'samples.tsv',tmp_path/'positions.tsv',tmp_path/'out',
            bin_genes=2,min_cpm=0,n_programmes=0,bin_stat=mode,counts_kind='raw',genome_build='hg38',
            coordinate_system='zero-based-half-open',gene_namespace='custom')
    assert pd.read_csv(out/'samples.tsv',sep='\t',dtype={'sample':str})['sample'].tolist() == names
    for field, name, dtype in ((0,'rna_state_posterior.npy',np.float32),(1,'rna_log_effect.npy',np.float32),
                               (2,'rna_state_index.npy',np.int8)):
        observed = np.load(out/name)
        assert observed.dtype == dtype
        assert np.array_equal(observed,np.stack([values[field] for values in captured]))
    segments = pd.read_csv(out/'segments.tsv',sep='\t')
    assert segments.n_bins.eq(3).all()
    assert not np.array_equal(np.load(out/'rna_log_effect.npy')[0],np.repeat(segments.rna_log_effect.iloc[0],3))
    metadata = json.loads((out/'run.json').read_text())['bin_rna_state_outputs']
    assert metadata['posterior']['axes'] == ['samples.tsv','bins.tsv','state_index']
    assert metadata['state_labels'] == ['loss','neutral','gain','amplification']
    assert metadata['effect']['units'] == 'natural_log_relative_rna_rate'
    assert metadata['additional_segment_averaging'] is False
    assert metadata['calibration_status'] == 'unvalidated' and metadata['biological_DNA_state_probability'] is False
    manifest = json.loads((out/'run_manifest.json').read_text())
    for name in ('rna_state_posterior.npy','rna_log_effect.npy','rna_state_index.npy'):
        assert name in manifest['output']['files']
