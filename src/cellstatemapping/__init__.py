"""
cellstatemapping — Reveal cancer-associated immune cell states.

A Python package for unsupervised machine learning and statistical modeling
of single-cell immune cell states, derived from Bachireddy et al. (2021)
Cell Reports.

Modules
-------
preprocessing
    Gene set and sample parsing utilities.
cfa
    Common Factor Analysis with Maximum Likelihood Estimation.
stats
    Weighted t-tests and significance evaluation.
distance
    Bhattacharyya distance and meta-clustering.
pipeline
    Automated pipeline connecting all modules.
cli
    Command-line interface.
"""

from cellstatemapping.preprocessing import (
    load_gene_set,
    load_curated_gene_sets,
    load_expression_matrix,
    normalize_log,
)
from cellstatemapping.cfa import (
    run_cfa,
    correlate_factors_with_geneset,
    permutation_test,
)
from cellstatemapping.stats import (
    weighted_t_test,
    weighted_cluster_comparison,
)
from cellstatemapping.distance import (
    bhattacharyya_distance,
    compute_distance_matrix,
)
from cellstatemapping.pipeline import CellStateMappingPipeline

__version__ = "1.0.0"

__all__ = [
    "load_gene_set",
    "load_curated_gene_sets",
    "load_expression_matrix",
    "normalize_log",
    "run_cfa",
    "correlate_factors_with_geneset",
    "permutation_test",
    "weighted_t_test",
    "weighted_cluster_comparison",
    "bhattacharyya_distance",
    "compute_distance_matrix",
    "CellStateMappingPipeline",
]
