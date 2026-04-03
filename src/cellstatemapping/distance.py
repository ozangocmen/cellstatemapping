"""
Distance metrics module for cellstatemapping.

Computes Bhattacharyya distance and pairwise distance matrices for
meta-clustering of cell clusters, matching the behavior of R's
``philentropy`` package.

Translates the Bhattacharyya distance computation from
``Objective 3.Rmd``.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def bhattacharyya_distance(
    p: np.ndarray,
    q: np.ndarray,
    epsilon: float = 1e-10,
) -> float:
    """Compute Bhattacharyya distance between two probability distributions.

    Matches the behavior and zero-handling of R's ``philentropy`` package,
    including epsilon padding for numerical stability.

    Parameters
    ----------
    p : np.ndarray
        First probability distribution.
    q : np.ndarray
        Second probability distribution.
    epsilon : float, default 1e-10
        Small constant to prevent log(0) and handle zero probabilities.

    Returns
    -------
    float
        Bhattacharyya distance.

    Notes
    -----
    The Bhattacharyya distance is defined as:
        D_B(p, q) = -ln(BC(p, q))
    where BC is the Bhattacharyya coefficient:
        BC(p, q) = sum(sqrt(p_i * q_i))

    R's philentropy uses natural log and epsilon padding, which this
    implementation replicates exactly.

    Examples
    --------
    >>> p = np.array([0.3, 0.3, 0.4])
    >>> q = np.array([0.2, 0.5, 0.3])
    >>> dist = bhattacharyya_distance(p, q)
    >>> print(f"{dist:.6f}")
    """
    p = np.asarray(p, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)

    if len(p) != len(q):
        raise ValueError(
            f"Distribution vectors must have equal length: {len(p)} != {len(q)}."
        )

    # Normalize to probability distributions
    p_sum = np.sum(p)
    q_sum = np.sum(q)

    if p_sum == 0 or q_sum == 0:
        raise ValueError("Cannot compute distance for zero-sum distributions.")

    p_norm = p / p_sum
    q_norm = q / q_sum

    # Bhattacharyya coefficient
    bc = np.sum(np.sqrt(p_norm * q_norm))

    # Clamp for numerical stability (matching philentropy behavior)
    bc = max(min(bc, 1.0), epsilon)

    # Bhattacharyya distance using natural log (philentropy default)
    return -np.log(bc)


def compute_distance_matrix(
    factor_matrix: pd.DataFrame,
) -> pd.DataFrame:
    """Compute pairwise Bhattacharyya distance matrix for meta-clustering.

    Parameters
    ----------
    factor_matrix : pd.DataFrame
        Matrix where columns represent clusters and rows represent factor
        loadings or features. Each column is treated as a distribution.

    Returns
    -------
    pd.DataFrame
        Symmetric distance matrix (clusters x clusters).

    Examples
    --------
    >>> dist_matrix = compute_distance_matrix(factor_loadings)
    >>> print(dist_matrix)
    """
    cols = factor_matrix.columns
    n = len(cols)
    dist_mat = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            # Ensure non-negative values for distance computation
            p = np.abs(factor_matrix.iloc[:, i].values)
            q = np.abs(factor_matrix.iloc[:, j].values)

            dist = bhattacharyya_distance(p, q)
            dist_mat[i, j] = dist
            dist_mat[j, i] = dist

    result = pd.DataFrame(dist_mat, index=cols, columns=cols)
    logger.info("Computed %d x %d pairwise distance matrix.", n, n)
    return result


def plot_distance_heatmap(
    distance_matrix: pd.DataFrame,
    output_path: Optional[str] = None,
    title: str = "Bhattacharyya Distance Heatmap",
    cmap: str = "viridis",
    figsize: tuple = (10, 8),
) -> None:
    """Plot a heatmap of the pairwise distance matrix.

    Parameters
    ----------
    distance_matrix : pd.DataFrame
        Symmetric distance matrix from ``compute_distance_matrix``.
    output_path : str or None
        If provided, saves the figure to this path.
    title : str
        Title for the plot.
    cmap : str
        Colormap for the heatmap.
    figsize : tuple
        Figure size in inches.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        distance_matrix,
        annot=True,
        fmt=".3f",
        cmap=cmap,
        square=True,
        ax=ax,
    )
    ax.set_title(title)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info("Saved distance heatmap to %s.", output_path)
    else:
        plt.show()

    plt.close(fig)
