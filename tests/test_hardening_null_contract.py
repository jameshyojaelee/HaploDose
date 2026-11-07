"""Unsupported observed/null mechanisms must fail before fitting."""
import pytest
from haplodose.clones import CloneParams, fit_specimen


def test_gene_veto_null_cannot_silently_omit_observed_evidence():
    def loader(_):
        pytest.fail('An unsupported experiment must not start reading gene evidence')
    with pytest.raises(ValueError, match='Gene-veto-enabled split null is unsupported'):
        fit_specimen(None, None, None, None, None, CloneParams(split_null_draws=19), gene_loader=loader)
