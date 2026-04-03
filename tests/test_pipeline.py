import pytest
from pathlib import Path

from cellstatemapping.pipeline import CellStateMappingPipeline

DATA_DIR = Path(__file__).parent.parent / "src" / "cellstatemapping" / "data"


def test_pipeline_end_to_end(tmp_path):
    output_dir = tmp_path / "results"
    
    pipeline = CellStateMappingPipeline(
        output_dir=str(output_dir),
        n_factors=2,
        n_permutations=10,
        n_bootstrap=0,
        seed=42
    )

    # 1. Load data
    pipeline.load_data(
        expression_files=[str(DATA_DIR / "test_expression_matrix.csv")],
        gene_set_file=str(DATA_DIR / "test_gene_set.txt"),
        metadata_file=str(DATA_DIR / "test_metadata.csv")
    )
    
    assert pipeline._gem is not None
    assert pipeline._gene_set is not None
    assert pipeline._metadata is not None

    # 2. Preprocess (no subsampling for small data)
    pipeline.preprocess(n_variable_features=50) # use 50 features
    assert pipeline._gem.shape[1] == 50

    # 3. Factor analysis
    pipeline.run_factor_analysis()
    assert pipeline._loadings is not None
    assert pipeline._correlations is not None

    # 4. Permutation test
    pipeline.run_permutation_tests()
    assert pipeline._perm_pvalues is not None
    assert "Factor1" in pipeline._perm_pvalues

    # 5. Weighted tests
    pipeline.run_weighted_tests()
    assert pipeline._test_results is not None

    # 6. Meta-clustering
    pipeline.compute_metaclusters()
    assert pipeline._distance_matrix is not None
    
    # 7. Save results
    pipeline.save_results()
    assert (output_dir / "factor_loadings.csv").exists()
    assert (output_dir / "factor_correlations.csv").exists()
    assert (output_dir / "permutation_pvalues.csv").exists()
    assert (output_dir / "weighted_test_results.csv").exists()
    assert (output_dir / "distance_matrix.csv").exists()
    assert (output_dir / "distance_heatmap.png").exists()
