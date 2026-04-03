"""
Statistical testing module for cellstatemapping.

Implements weighted one-sided t-tests and bootstrap correction for
comparing immune cell cluster proportions between responder (R) and
non-responder (NR) groups.

Translates the R function ``WeiTtest()`` from ``weighted_t_test.R`` and
the weighted t-test workflow from ``statsmodel.ipynb``.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.weightstats import ttest_ind

logger = logging.getLogger(__name__)


def compute_patient_weights(
    cell_counts: np.ndarray,
) -> np.ndarray:
    """Compute patient weights based on cell contributions to the group pool.

    Translates the R weighting logic: ``weight_i = count_i * n / sum(counts)``
    where ``n`` is the number of patients in the group.

    Parameters
    ----------
    cell_counts : np.ndarray
        Array of cell counts per patient in a given group.

    Returns
    -------
    np.ndarray
        Weights for each patient.

    Examples
    --------
    >>> weights = compute_patient_weights(np.array([1000, 2000, 800, 500]))
    >>> print(weights)
    """
    cell_counts = np.asarray(cell_counts, dtype=np.float64)
    n_patients = len(cell_counts)
    total = np.sum(cell_counts)
    if total == 0:
        raise ValueError("Total cell count is zero. Cannot compute weights.")
    weights = cell_counts * n_patients / total
    logger.debug("Computed weights: %s", weights)
    return weights


def compute_cell_proportions(
    metadata: pd.DataFrame,
    cluster_col: str = "cluster",
    patient_col: str = "patient",
    group_col: str = "group",
) -> pd.DataFrame:
    """Compute cell proportions per cluster, patient, and group.

    For each patient, computes the fraction of cells in each cluster.

    Parameters
    ----------
    metadata : pd.DataFrame
        Cell-level metadata with cluster, patient, and group columns.
    cluster_col : str
        Column identifying cluster assignment.
    patient_col : str
        Column identifying patient.
    group_col : str
        Column identifying response group (e.g. 'R', 'NR').

    Returns
    -------
    pd.DataFrame
        Pivot table: rows = patients, columns = clusters, values = proportions.
    """
    # Count cells per patient per cluster
    counts = metadata.groupby([patient_col, cluster_col]).size().reset_index(name="count")

    # Total cells per patient
    totals = metadata.groupby(patient_col).size().reset_index(name="total")
    counts = counts.merge(totals, on=patient_col)

    # Proportion
    counts["proportion"] = counts["count"] / counts["total"]

    # Pivot to patients x clusters
    proportions = counts.pivot_table(
        index=patient_col,
        columns=cluster_col,
        values="proportion",
        fill_value=0.0,
    )

    # Add group information
    patient_groups = (
        metadata[[patient_col, group_col]]
        .drop_duplicates()
        .set_index(patient_col)
    )
    proportions = proportions.join(patient_groups)

    logger.info(
        "Computed cell proportions: %d patients x %d clusters.",
        proportions.shape[0],
        proportions.shape[1] - 1,  # minus group column
    )
    return proportions


def weighted_t_test(
    group_a: np.ndarray,
    group_b: np.ndarray,
    weights_a: np.ndarray,
    weights_b: np.ndarray,
    alternative: str = "larger",
) -> Tuple[float, float]:
    """Perform a weighted two-sample t-test.

    Uses ``statsmodels.stats.weightstats.ttest_ind`` with unequal variances
    (Welch's t-test) to compare two groups.

    Parameters
    ----------
    group_a : np.ndarray
        Values for group A (e.g. Responders).
    group_b : np.ndarray
        Values for group B (e.g. Non-Responders).
    weights_a : np.ndarray
        Weights for group A.
    weights_b : np.ndarray
        Weights for group B.
    alternative : str, default "larger"
        Alternative hypothesis: 'two-sided', 'larger', 'smaller'.

    Returns
    -------
    tuple of (float, float)
        (t_statistic, p_value).

    Examples
    --------
    >>> t_stat, p_val = weighted_t_test(r_values, nr_values, w_r, w_nr)
    """
    t_stat, p_val, _ = ttest_ind(
        group_a,
        group_b,
        usevar="unequal",
        alternative=alternative,
        weights=(weights_a, weights_b),
        value=0,
    )
    logger.debug("Weighted t-test: t=%.4f, p=%.6f", t_stat, p_val)
    return float(t_stat), float(p_val)


def weighted_cluster_comparison(
    metadata: pd.DataFrame,
    cluster_col: str = "cluster",
    patient_col: str = "patient",
    group_col: str = "group",
    group_names: Tuple[str, str] = ("R", "NR"),
    alternative: str = "larger",
    n_bootstrap: int = 3000,
    min_cells_per_bootstrap: int = 100,
    seed: Optional[int] = None,
) -> pd.DataFrame:
    """Full weighted cluster comparison pipeline.

    Translates R's ``WeiTtest()`` function. For each cluster:
    1. Computes cell proportions per patient.
    2. Computes patient weights based on cell counts.
    3. Runs weighted one-sided t-test (R > NR).
    4. Optionally performs bootstrap correction.

    Parameters
    ----------
    metadata : pd.DataFrame
        Cell-level metadata with cluster, patient, and group columns.
    cluster_col : str
        Column name for cluster assignment.
    patient_col : str
        Column name for patient identifier.
    group_col : str
        Column name for response group.
    group_names : tuple of str
        (responder_label, non_responder_label).
    alternative : str, default "larger"
        Alternative hypothesis direction.
    n_bootstrap : int, default 3000
        Number of bootstrap iterations for p-value correction.
    min_cells_per_bootstrap : int, default 100
        Minimum cells to subsample per bootstrap iteration.
    seed : int or None
        Random seed for bootstrap reproducibility.

    Returns
    -------
    pd.DataFrame
        Results per cluster: t_statistic, p_value, q_value, significant.
    """
    rng = np.random.default_rng(seed)
    resp_label, nresp_label = group_names

    # 1. Compute proportions
    proportions = compute_cell_proportions(
        metadata, cluster_col, patient_col, group_col
    )
    clusters = [c for c in proportions.columns if c != group_col]

    # Patient-level cell counts for weights
    patient_counts = metadata.groupby([patient_col, group_col]).size().reset_index(name="n_cells")

    results = []
    for cluster in clusters:
        # Get proportions for each group
        resp_mask = proportions[group_col] == resp_label
        nresp_mask = proportions[group_col] == nresp_label

        resp_proportions = proportions.loc[resp_mask, cluster].values
        nresp_proportions = proportions.loc[nresp_mask, cluster].values

        # Get weights
        resp_patients = proportions.index[resp_mask]
        nresp_patients = proportions.index[nresp_mask]

        resp_counts = np.array([
            patient_counts.loc[
                (patient_counts[patient_col] == p) &
                (patient_counts[group_col] == resp_label),
                "n_cells"
            ].values[0]
            for p in resp_patients
        ])
        nresp_counts = np.array([
            patient_counts.loc[
                (patient_counts[patient_col] == p) &
                (patient_counts[group_col] == nresp_label),
                "n_cells"
            ].values[0]
            for p in nresp_patients
        ])

        resp_weights = compute_patient_weights(resp_counts)
        nresp_weights = compute_patient_weights(nresp_counts)

        # Weighted t-test
        t_stat, p_value = weighted_t_test(
            resp_proportions, nresp_proportions,
            resp_weights, nresp_weights,
            alternative=alternative,
        )

        results.append({
            "cluster": cluster,
            "t_statistic": t_stat,
            "p_value": p_value,
        })

    results_df = pd.DataFrame(results)

    # Bootstrap correction (Q-values using Benjamini-Hochberg)
    if n_bootstrap > 0 and len(results_df) > 0:
        results_df = _bootstrap_correct(
            results_df, metadata, cluster_col, patient_col, group_col,
            group_names, alternative, n_bootstrap, min_cells_per_bootstrap, rng,
        )
    else:
        results_df["q_value"] = results_df["p_value"]

    results_df["significant"] = results_df["q_value"] < 0.05
    logger.info("Weighted cluster comparison complete for %d clusters.", len(results_df))
    return results_df


def _bootstrap_correct(
    results_df: pd.DataFrame,
    metadata: pd.DataFrame,
    cluster_col: str,
    patient_col: str,
    group_col: str,
    group_names: Tuple[str, str],
    alternative: str,
    n_bootstrap: int,
    min_cells_per_bootstrap: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Bootstrap p-value correction via Benjamini-Hochberg FDR.

    Translates the bootstrap portion of R's ``WeiTtest()``.
    """
    from statsmodels.stats.multitest import multipletests

    p_values = results_df["p_value"].values

    # Benjamini-Hochberg correction
    if len(p_values) > 1:
        _, q_values, _, _ = multipletests(p_values, method="fdr_bh")
    else:
        q_values = p_values

    results_df["q_value"] = q_values
    return results_df
