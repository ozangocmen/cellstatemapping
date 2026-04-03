"""
Common Factor Analysis (CFA) module for cellstatemapping.

Replicates R's ``factanal()`` with Maximum Likelihood Estimation (MLE)
and varimax rotation using the ``factor_analyzer`` package. Provides
permutation testing for factor–gene-set correlation significance.

This module translates the following R functions:
- ``RunCFA_GScors()`` from ``BIGQA-WP1-functions.R``
- ``permutationTest()`` from ``BIGQA-WP1-functions.R``
"""

from __future__ import annotations

import logging
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from factor_analyzer import FactorAnalyzer
from scipy.stats import pearsonr

logger = logging.getLogger(__name__)


def run_cfa(
    gem: pd.DataFrame,
    n_factors: int = 3,
    rotation: str = "varimax",
    method: str = "ml",
) -> pd.DataFrame:
    """Perform Common Factor Analysis with MLE and varimax rotation.

    Replicates R's ``factanal(GEM, nF, rotation="varimax")`` using the
    ``factor_analyzer`` package with Maximum Likelihood ('ml') estimation.

    Parameters
    ----------
    gem : pd.DataFrame
        Gene expression matrix with shape (n_cells, n_genes).
        Rows are observations (cells), columns are variables (genes).
    n_factors : int, default 3
        Number of factors to extract.
    rotation : str, default "varimax"
        Rotation method (e.g. 'varimax', 'promax').
    method : str, default "ml"
        Extraction method. 'ml' for Maximum Likelihood (matches R's factanal).

    Returns
    -------
    pd.DataFrame
        Factor loadings matrix with shape (n_genes, n_factors).
        Index = gene names, Columns = Factor1, Factor2, ...

    Examples
    --------
    >>> loadings = run_cfa(gem, n_factors=3)
    >>> print(loadings.head())
    """
    fa = FactorAnalyzer(n_factors=n_factors, rotation=rotation, method=method)
    fa.fit(gem)

    factor_names = [f"Factor{i + 1}" for i in range(n_factors)]
    loadings = pd.DataFrame(
        fa.loadings_,
        index=gem.columns,
        columns=factor_names,
    )
    logger.info(
        "CFA complete: %d factors extracted from %d genes x %d cells.",
        n_factors,
        gem.shape[1],
        gem.shape[0],
    )
    return loadings


def _compute_mean_geneset_expression(
    gem: pd.DataFrame,
    gene_set: list,
) -> pd.Series:
    """Compute the mean expression of a gene set per cell.

    Parameters
    ----------
    gem : pd.DataFrame
        Gene expression matrix (cells x genes).
    gene_set : list of str
        Gene symbols to average.

    Returns
    -------
    pd.Series
        Mean gene set expression per cell.
    """
    overlap = [g for g in gene_set if g in gem.columns]
    if not overlap:
        raise ValueError("No genes from the gene set are found in the expression matrix.")
    logger.info(
        "Gene set overlap: %d / %d genes present in the expression matrix.",
        len(overlap),
        len(gene_set),
    )
    return gem[overlap].mean(axis=1)


def correlate_factors_with_geneset(
    gem: pd.DataFrame,
    gene_set: list,
    n_factors: int = 3,
    rotation: str = "varimax",
    method: str = "ml",
) -> Tuple[pd.DataFrame, Dict[str, Tuple[float, float]]]:
    """Run CFA and correlate each factor with mean gene set expression.

    This is the Python equivalent of R's ``RunCFA_GScors()`` function.

    Parameters
    ----------
    gem : pd.DataFrame
        Gene expression matrix (cells x genes).
    gene_set : list of str
        Gene symbols for correlation.
    n_factors : int, default 3
        Number of factors to extract.
    rotation : str, default "varimax"
        Rotation method.
    method : str, default "ml"
        Factor extraction method.

    Returns
    -------
    loadings : pd.DataFrame
        Factor loadings matrix (genes x factors).
    correlations : dict
        Dictionary mapping factor name -> (correlation, p_value).

    Examples
    --------
    >>> loadings, corrs = correlate_factors_with_geneset(gem, gene_set)
    >>> for factor, (r, p) in corrs.items():
    ...     print(f"{factor}: r={r:.4f}, p={p:.6f}")
    """
    loadings = run_cfa(gem, n_factors=n_factors, rotation=rotation, method=method)
    mean_gs_expr = _compute_mean_geneset_expression(gem, gene_set)

    # Factor scores for correlation: cells x factors
    fa = FactorAnalyzer(n_factors=n_factors, rotation=rotation, method=method)
    fa.fit(gem)
    factor_scores = pd.DataFrame(
        fa.transform(gem),
        index=gem.index,
        columns=loadings.columns,
    )

    correlations: Dict[str, Tuple[float, float]] = {}
    for col in factor_scores.columns:
        r, p = pearsonr(factor_scores[col], mean_gs_expr)
        correlations[col] = (r, p)
        logger.info("Correlation %s with gene set: r=%.4f, p=%.6f", col, r, p)

    return loadings, correlations


def permutation_test(
    observed_corr: float,
    gem: pd.DataFrame,
    gene_set: list,
    factor_scores: pd.Series,
    n_permutations: int = 500,
    seed: int = 878712,
) -> float:
    """Permutation test for factor–gene-set correlation significance.

    Translates R's ``permutationTest()`` from ``BIGQA-WP1-functions.R``.
    Randomly selects a gene set of the same size from all genes, computes
    the correlation with factor scores, and builds a null distribution.

    Parameters
    ----------
    observed_corr : float
        Observed Pearson correlation between factor scores and
        mean gene set expression.
    gem : pd.DataFrame
        Gene expression matrix (cells x genes).
    gene_set : list of str
        Original gene set used to compute observed_corr.
    factor_scores : pd.Series
        Factor scores for a single factor (indexed by cell).
    n_permutations : int, default 500
        Number of random gene sets to draw.
    seed : int, default 878712
        Random seed for reproducibility (matches R default).

    Returns
    -------
    float
        Empirical p-value: fraction of permuted correlations >= observed.

    Notes
    -----
    The R implementation draws gene sets of the same size as ``gene_set``
    from all genes in ``gem``, computes mean expression per cell, then
    correlates that with ``factor_scores``.
    """
    rng = np.random.default_rng(seed)
    all_genes = list(gem.columns)
    n_genes = len([g for g in gene_set if g in all_genes])

    if n_genes == 0:
        raise ValueError("No overlapping genes between gene set and expression matrix.")

    null_corrs = np.zeros(n_permutations)
    for i in range(n_permutations):
        random_genes = rng.choice(all_genes, size=n_genes, replace=False)
        random_mean = gem[random_genes].mean(axis=1)
        r, _ = pearsonr(factor_scores, random_mean)
        null_corrs[i] = abs(r)

    p_value = np.mean(null_corrs >= abs(observed_corr))
    logger.info(
        "Permutation test: observed_corr=%.4f, p_value=%.6f (%d permutations)",
        observed_corr,
        p_value,
        n_permutations,
    )
    return p_value
