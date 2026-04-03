import pytest
import pandas as pd
import numpy as np

from cellstatemapping.cfa import (
    run_cfa,
    _compute_mean_geneset_expression,
    correlate_factors_with_geneset,
    permutation_test
)


@pytest.fixture
def sample_gem():
    # 50 cells, 20 genes
    np.random.seed(42)
    data = np.random.rand(50, 20)
    genes = [f"Gene{i}" for i in range(1, 21)]
    cells = [f"Cell{i}" for i in range(1, 51)]
    return pd.DataFrame(data, index=cells, columns=genes)


@pytest.fixture
def sample_gene_set():
    return ["Gene1", "Gene2", "Gene3"]


def test_run_cfa(sample_gem):
    loadings = run_cfa(sample_gem, n_factors=2, method="ml")
    assert loadings.shape == (20, 2)
    assert "Factor1" in loadings.columns
    assert "Factor2" in loadings.columns


def test_compute_mean_geneset_expression(sample_gem, sample_gene_set):
    mean_expr = _compute_mean_geneset_expression(sample_gem, sample_gene_set)
    assert len(mean_expr) == 50
    # Manually calculate for Cell1
    val = sample_gem.loc["Cell1", ["Gene1", "Gene2", "Gene3"]].mean()
    assert np.isclose(mean_expr.loc["Cell1"], val)


def test_correlate_factors_with_geneset(sample_gem, sample_gene_set):
    loadings, corrs = correlate_factors_with_geneset(
        sample_gem, sample_gene_set, n_factors=2
    )
    assert "Factor1" in corrs
    assert "Factor2" in corrs
    # (correlation, p_value)
    assert len(corrs["Factor1"]) == 2


def test_permutation_test(sample_gem, sample_gene_set):
    # Dummy factor scores (50 cells)
    np.random.seed(42)
    factor_scores = pd.Series(np.random.rand(50), index=sample_gem.index)
    
    # Calculate observed corr
    mean_expr = _compute_mean_geneset_expression(sample_gem, sample_gene_set)
    observed_corr = np.corrcoef(factor_scores, mean_expr)[0, 1]
    
    p_val = permutation_test(
        observed_corr=observed_corr,
        gem=sample_gem,
        gene_set=sample_gene_set,
        factor_scores=factor_scores,
        n_permutations=20,
        seed=42
    )
    assert 0.0 <= p_val <= 1.0
