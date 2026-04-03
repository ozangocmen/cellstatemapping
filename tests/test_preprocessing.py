import pytest
from pathlib import Path
import pandas as pd
import numpy as np

from cellstatemapping.preprocessing import (
    load_gene_set,
    load_curated_gene_sets,
    load_expression_matrix,
    normalize_log,
    find_variable_features,
    subsample_cells
)

DATA_DIR = Path(__file__).parent.parent / "src" / "cellstatemapping" / "data"


def test_load_gene_set():
    gene_set_path = DATA_DIR / "test_gene_set.txt"
    genes = load_gene_set(gene_set_path)
    assert len(genes) == 20
    assert "Gene1" in genes


def test_load_curated_gene_sets():
    curated_path = DATA_DIR / "test_curated_gene_sets.txt"
    gene_sets = load_curated_gene_sets(curated_path)
    assert "GeneSet1" in gene_sets
    assert "GeneSet2" in gene_sets
    assert len(gene_sets["GeneSet1"]) == 10
    assert len(gene_sets["GeneSet2"]) == 10


def test_load_expression_matrix():
    gem_path = DATA_DIR / "test_expression_matrix.csv"
    gem = load_expression_matrix(gem_path, transpose=True)
    # Original is 200 cells x 100 genes. Transposed: 100 genes x 200 cells.
    assert gem.shape == (100, 200)
    assert "Cell1" in gem.columns
    assert "Gene1" in gem.index


def test_subsample_cells():
    gem_path = DATA_DIR / "test_expression_matrix.csv"
    gem = load_expression_matrix(gem_path, transpose=True)
    subsampled = subsample_cells(gem, n_cells=50, random_state=42)
    assert subsampled.shape == (100, 50)
    
    # Subsampling more cells than available should return all
    all_cells = subsample_cells(gem, n_cells=500)
    assert all_cells.shape == (100, 200)


def test_normalize_log():
    # Small test matrix
    df_raw = pd.DataFrame({
        "cell1": [10, 20, 0],
        "cell2": [0, 50, 50]
    }, index=["gene1", "gene2", "gene3"])
    
    df_norm = normalize_log(df_raw, scale_factor=100)
    
    # Check cell1: total=30. gene1 = ln(1 + 10/30 * 100) = ln(1 + 33.33)
    val = np.log1p(10/30 * 100)
    assert np.isclose(df_norm.loc["gene1", "cell1"], val)


def test_find_variable_features():
    df = pd.DataFrame({
        "c1": [10, 10, 1],
        "c2": [10, 20, 1],
        "c3": [10, 30, 2],
        "c4": [10, 40, 1]
    }, index=["g1", "g2", "g3"])
    
    # g2 has highest variance
    var_genes = find_variable_features(df, n_features=2)
    assert "g2" in var_genes
    assert len(var_genes) == 2
