"""
Automated pipeline for cellstatemapping.

Connects preprocessing, common factor analysis, statistical testing,
and meta-clustering into a single end-to-end workflow. Translates
``final_script.Rmd`` into a reusable Python pipeline class.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from cellstatemapping.preprocessing import (
    load_expression_matrix,
    load_gene_set,
    load_curated_gene_sets,
    normalize_log,
    find_variable_features,
    subsample_cells,
)
from cellstatemapping.cfa import (
    run_cfa,
    correlate_factors_with_geneset,
    permutation_test,
)
from cellstatemapping.stats import (
    weighted_cluster_comparison,
)
from cellstatemapping.distance import (
    compute_distance_matrix,
    plot_distance_heatmap,
)

logger = logging.getLogger(__name__)


class CellStateMappingPipeline:
    """End-to-end pipeline for cancer-associated immune cell state mapping.

    This class orchestrates the full analytical workflow:

    1. **Load data**: Expression matrices and gene sets
    2. **Preprocess**: Log-normalization, variable feature selection
    3. **Factor analysis**: CFA with MLE + varimax, gene set correlation
    4. **Permutation tests**: Significance of factor–gene-set associations
    5. **Weighted t-tests**: Responder vs Non-Responder cluster comparisons
    6. **Meta-clustering**: Bhattacharyya distance-based cluster similarity

    Parameters
    ----------
    output_dir : str or Path
        Directory for saving results and figures.
    n_factors : int, default 3
        Number of factors for CFA.
    n_permutations : int, default 500
        Number of permutations for significance testing.
    n_bootstrap : int, default 3000
        Number of bootstrap iterations for p-value correction.
    scale_factor : float, default 1e4
        Scale factor for log-normalization.
    seed : int, default 42
        Random seed for reproducibility.

    Examples
    --------
    >>> pipeline = CellStateMappingPipeline(output_dir="results/")
    >>> pipeline.load_data(
    ...     expression_files=["sample1.csv", "sample2.csv"],
    ...     gene_set_file="GS.proinf.txt",
    ...     metadata_file="metadata.csv",
    ... )
    >>> pipeline.preprocess()
    >>> pipeline.run_factor_analysis()
    >>> pipeline.run_permutation_tests()
    >>> pipeline.run_weighted_tests(group_names=("R", "NR"))
    >>> results = pipeline.get_results()
    """

    def __init__(
        self,
        output_dir: str = "results",
        n_factors: int = 3,
        n_permutations: int = 500,
        n_bootstrap: int = 3000,
        scale_factor: float = 1e4,
        seed: int = 42,
    ) -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.n_factors = n_factors
        self.n_permutations = n_permutations
        self.n_bootstrap = n_bootstrap
        self.scale_factor = scale_factor
        self.seed = seed

        # Internal state
        self._gem: Optional[pd.DataFrame] = None
        self._gene_set: Optional[List[str]] = None
        self._metadata: Optional[pd.DataFrame] = None
        self._loadings: Optional[pd.DataFrame] = None
        self._correlations: Optional[Dict[str, Tuple[float, float]]] = None
        self._perm_pvalues: Optional[Dict[str, float]] = None
        self._test_results: Optional[pd.DataFrame] = None
        self._distance_matrix: Optional[pd.DataFrame] = None

        logger.info("Pipeline initialized (output_dir=%s).", self.output_dir)

    def load_data(
        self,
        expression_files: Optional[List[str]] = None,
        gene_set_file: Optional[str] = None,
        curated_gene_set_file: Optional[str] = None,
        metadata_file: Optional[str] = None,
        gem: Optional[pd.DataFrame] = None,
        metadata: Optional[pd.DataFrame] = None,
        gene_set: Optional[List[str]] = None,
    ) -> "CellStateMappingPipeline":
        """Load expression matrices, gene sets, and metadata.

        Can accept file paths or pre-loaded DataFrames.

        Parameters
        ----------
        expression_files : list of str, optional
            Paths to CSV expression matrices.
        gene_set_file : str, optional
            Path to single gene set file.
        curated_gene_set_file : str, optional
            Path to curated (multi) gene set file.
        metadata_file : str, optional
            Path to metadata CSV.
        gem : pd.DataFrame, optional
            Pre-loaded expression matrix (cells x genes).
        metadata : pd.DataFrame, optional
            Pre-loaded metadata DataFrame.
        gene_set : list of str, optional
            Pre-loaded gene set.

        Returns
        -------
        self
        """
        # Load gene set
        if gene_set is not None:
            self._gene_set = gene_set
        elif gene_set_file is not None:
            self._gene_set = load_gene_set(gene_set_file)
        elif curated_gene_set_file is not None:
            gene_sets = load_curated_gene_sets(curated_gene_set_file)
            # Use the first gene set by default
            first_key = next(iter(gene_sets))
            self._gene_set = gene_sets[first_key]
            logger.info("Using gene set '%s' from curated file.", first_key)

        # Load expression matrix
        if gem is not None:
            self._gem = gem
        elif expression_files is not None:
            dfs = []
            for fp in expression_files:
                df = load_expression_matrix(fp, transpose=True)
                dfs.append(df)
            if len(dfs) == 1:
                self._gem = dfs[0].T  # back to cells x genes
            else:
                # Concatenate along cell axis, align genes
                combined = pd.concat([d.T for d in dfs], axis=0, join="inner")
                self._gem = combined
            logger.info("Loaded GEM: %d cells x %d genes.", *self._gem.shape)

        # Load metadata
        if metadata is not None:
            self._metadata = metadata
        elif metadata_file is not None:
            self._metadata = pd.read_csv(metadata_file, index_col=0)

        return self

    def preprocess(
        self,
        n_variable_features: int = 2000,
        subsample_n: Optional[int] = None,
    ) -> "CellStateMappingPipeline":
        """Preprocess the expression matrix.

        Performs log-normalization and optionally selects variable features
        and subsamples cells.

        Parameters
        ----------
        n_variable_features : int, default 2000
            Number of variable features to select.
        subsample_n : int or None
            If provided, subsample this many cells.

        Returns
        -------
        self
        """
        if self._gem is None:
            raise ValueError("No expression matrix loaded. Call load_data() first.")

        # Expression matrix is cells x genes
        gem_t = self._gem.T  # genes x cells

        # Log-normalize
        gem_norm = normalize_log(gem_t, scale_factor=self.scale_factor)

        # Variable features
        var_features = find_variable_features(gem_norm, n_features=n_variable_features)
        gem_norm = gem_norm.loc[var_features]

        # Back to cells x genes
        self._gem = gem_norm.T

        # Subsample
        if subsample_n is not None:
            self._gem = subsample_cells(
                self._gem.T, n_cells=subsample_n, random_state=self.seed
            ).T

        logger.info("Preprocessing complete: %d cells x %d genes.", *self._gem.shape)
        return self

    def run_factor_analysis(self) -> "CellStateMappingPipeline":
        """Run CFA and correlate factors with gene set.

        Returns
        -------
        self
        """
        if self._gem is None:
            raise ValueError("No expression matrix. Call load_data()/preprocess().")
        if self._gene_set is None:
            raise ValueError("No gene set loaded. Call load_data().")

        self._loadings, self._correlations = correlate_factors_with_geneset(
            self._gem, self._gene_set, n_factors=self.n_factors
        )
        logger.info("Factor analysis complete.")
        return self

    def run_permutation_tests(self) -> "CellStateMappingPipeline":
        """Run permutation tests for each factor–gene-set correlation.

        Returns
        -------
        self
        """
        if self._correlations is None or self._loadings is None:
            raise ValueError("Run factor analysis first.")

        from factor_analyzer import FactorAnalyzer

        fa = FactorAnalyzer(
            n_factors=self.n_factors, rotation="varimax", method="ml"
        )
        fa.fit(self._gem)
        factor_scores = pd.DataFrame(
            fa.transform(self._gem),
            index=self._gem.index,
            columns=self._loadings.columns,
        )

        self._perm_pvalues = {}
        for factor_name, (corr, _) in self._correlations.items():
            p = permutation_test(
                observed_corr=corr,
                gem=self._gem,
                gene_set=self._gene_set,
                factor_scores=factor_scores[factor_name],
                n_permutations=self.n_permutations,
                seed=self.seed,
            )
            self._perm_pvalues[factor_name] = p

        logger.info("Permutation tests complete.")
        return self

    def run_weighted_tests(
        self,
        group_names: Tuple[str, str] = ("R", "NR"),
        cluster_col: str = "cluster",
        patient_col: str = "patient",
        group_col: str = "group",
    ) -> "CellStateMappingPipeline":
        """Run weighted t-tests comparing groups across clusters.

        Parameters
        ----------
        group_names : tuple of str
            (responder_label, non_responder_label).
        cluster_col, patient_col, group_col : str
            Column names in metadata.

        Returns
        -------
        self
        """
        if self._metadata is None:
            raise ValueError("No metadata loaded. Call load_data().")

        self._test_results = weighted_cluster_comparison(
            metadata=self._metadata,
            cluster_col=cluster_col,
            patient_col=patient_col,
            group_col=group_col,
            group_names=group_names,
            n_bootstrap=self.n_bootstrap,
            seed=self.seed,
        )
        logger.info("Weighted tests complete.")
        return self

    def compute_metaclusters(self) -> "CellStateMappingPipeline":
        """Compute Bhattacharyya distance matrix for meta-clustering.

        Returns
        -------
        self
        """
        if self._loadings is None:
            raise ValueError("Run factor analysis first.")

        self._distance_matrix = compute_distance_matrix(self._loadings)

        # Save heatmap
        plot_distance_heatmap(
            self._distance_matrix,
            output_path=str(self.output_dir / "distance_heatmap.png"),
        )
        logger.info("Meta-clustering distance matrix computed.")
        return self

    def get_results(self) -> Dict[str, Any]:
        """Return all computed results.

        Returns
        -------
        dict
            Dictionary with keys: 'loadings', 'correlations',
            'permutation_pvalues', 'test_results', 'distance_matrix'.
        """
        return {
            "loadings": self._loadings,
            "correlations": self._correlations,
            "permutation_pvalues": self._perm_pvalues,
            "test_results": self._test_results,
            "distance_matrix": self._distance_matrix,
        }

    def save_results(self) -> None:
        """Save all computed results to the output directory."""
        if self._loadings is not None:
            self._loadings.to_csv(self.output_dir / "factor_loadings.csv")
        if self._test_results is not None:
            self._test_results.to_csv(
                self.output_dir / "weighted_test_results.csv", index=False
            )
        if self._distance_matrix is not None:
            self._distance_matrix.to_csv(self.output_dir / "distance_matrix.csv")
        if self._correlations is not None:
            corr_df = pd.DataFrame(
                {k: {"correlation": v[0], "p_value": v[1]}
                 for k, v in self._correlations.items()}
            ).T
            corr_df.to_csv(self.output_dir / "factor_correlations.csv")
        if self._perm_pvalues is not None:
            perm_df = pd.DataFrame(
                list(self._perm_pvalues.items()),
                columns=["factor", "permutation_p_value"],
            )
            perm_df.to_csv(
                self.output_dir / "permutation_pvalues.csv", index=False
            )
        logger.info("Results saved to %s.", self.output_dir)
