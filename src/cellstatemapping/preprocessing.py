"""
Preprocessing module for cellstatemapping.

Provides utilities for loading gene expression matrices, parsing gene sets,
normalizing data, and subsampling cells. These functions replicate the
preprocessing steps from Seurat/R workflows in pure Python.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def load_gene_set(filepath: Union[str, Path]) -> List[str]:
    """Load a gene set from a single-column text file.

    Each line in the file should contain one gene symbol.
    Blank lines and duplicate entries are automatically removed.

    Parameters
    ----------
    filepath : str or Path
        Path to a text file with one gene symbol per line.

    Returns
    -------
    list of str
        Sorted, deduplicated list of gene symbols.

    Examples
    --------
    >>> genes = load_gene_set("GS.proinf.txt")
    >>> print(genes[:3])
    ['CCL2', 'CCL20', 'CCL3']
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Gene set file not found: {filepath}")

    genes: List[str] = []
    with open(filepath, "r") as f:
        for line in f:
            gene = line.strip()
            if gene and gene not in genes:
                genes.append(gene)

    logger.info("Loaded %d genes from %s", len(genes), filepath.name)
    return sorted(genes)


def load_curated_gene_sets(filepath: Union[str, Path]) -> dict:
    """Load curated gene sets from a tab-delimited file.

    The file format should have each row as a gene set, with columns:
    - Column 0: Gene set ID
    - Column 1: Gene set description
    - Columns 2+: Gene symbols

    Parameters
    ----------
    filepath : str or Path
        Path to the tab-delimited gene set file.

    Returns
    -------
    dict
        Dictionary mapping gene set ID to list of gene symbols.

    Examples
    --------
    >>> gene_sets = load_curated_gene_sets("CuratedGeneSets.txt")
    >>> print(gene_sets.keys())
    dict_keys(['ProInf', 'CeP'])
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Curated gene set file not found: {filepath}")

    gene_sets: dict = {}
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                gs_id = parts[0]
                # parts[1] is the description, skip it
                genes = [g.strip() for g in parts[2:] if g.strip()]
                gene_sets[gs_id] = sorted(set(genes))
                logger.info(
                    "Loaded gene set '%s' with %d genes", gs_id, len(gene_sets[gs_id])
                )

    return gene_sets


def load_expression_matrix(
    filepath: Union[str, Path],
    transpose: bool = True,
) -> pd.DataFrame:
    """Load a gene expression matrix from a CSV file.

    By default, the input is expected to have cells as rows and genes as columns
    (the standard dense CSV format from Bachireddy et al.), and will be transposed
    so that genes are rows and cells are columns.

    Parameters
    ----------
    filepath : str or Path
        Path to the CSV file containing the expression matrix.
    transpose : bool, default True
        If True, transpose the matrix so genes are rows and cells are columns.

    Returns
    -------
    pd.DataFrame
        Gene expression matrix with genes as rows and cells as columns.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Expression matrix file not found: {filepath}")

    df = pd.read_csv(filepath, index_col=0, header=0)
    if transpose:
        df = df.T

    logger.info(
        "Loaded expression matrix: %d genes x %d cells from %s",
        df.shape[0],
        df.shape[1],
        filepath.name,
    )
    return df


def subsample_cells(
    gem: pd.DataFrame,
    n_cells: int = 1000,
    random_state: Optional[int] = None,
) -> pd.DataFrame:
    """Randomly subsample cells from a gene expression matrix.

    Parameters
    ----------
    gem : pd.DataFrame
        Gene expression matrix (genes x cells).
    n_cells : int, default 1000
        Number of cells to retain.
    random_state : int or None, default None
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        Subsampled gene expression matrix.
    """
    rng = np.random.default_rng(random_state)
    n_available = gem.shape[1]

    if n_cells >= n_available:
        logger.warning(
            "Requested %d cells but only %d available; returning all.",
            n_cells,
            n_available,
        )
        return gem

    selected = rng.choice(n_available, size=n_cells, replace=False)
    selected.sort()
    result = gem.iloc[:, selected]
    logger.info("Subsampled %d cells from %d total.", n_cells, n_available)
    return result


def normalize_log(
    gem: pd.DataFrame,
    scale_factor: float = 1e4,
) -> pd.DataFrame:
    """Log-normalize a gene expression matrix.

    Replicates Seurat's LogNormalize method:
    ``log1p(count / total_count_per_cell * scale_factor)``

    Parameters
    ----------
    gem : pd.DataFrame
        Raw count gene expression matrix (genes x cells).
    scale_factor : float, default 1e4
        Scale factor applied before log transformation.

    Returns
    -------
    pd.DataFrame
        Log-normalized gene expression matrix.
    """
    # Column sums = total counts per cell
    col_sums = gem.sum(axis=0)
    # Avoid division by zero
    col_sums = col_sums.replace(0, 1)

    normalized = gem.div(col_sums, axis=1) * scale_factor
    normalized = np.log1p(normalized)

    logger.info("Log-normalized expression matrix (scale_factor=%g).", scale_factor)
    return normalized


def find_variable_features(
    gem: pd.DataFrame,
    n_features: int = 2000,
) -> List[str]:
    """Select highly variable genes using variance-based selection.

    A simplified version of Seurat's FindVariableFeatures with 'vst' method.
    Selects the top ``n_features`` genes by variance across cells.

    Parameters
    ----------
    gem : pd.DataFrame
        Gene expression matrix (genes x cells).
    n_features : int, default 2000
        Number of top variable features to select.

    Returns
    -------
    list of str
        Gene names of the most variable features.
    """
    variances = gem.var(axis=1)
    top_genes = variances.nlargest(n_features).index.tolist()
    logger.info("Selected %d variable features.", len(top_genes))
    return top_genes
