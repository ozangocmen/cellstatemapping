# Cell State Mapp
<div align="center">
<img src="logo.png" style="background-color: white; padding: 10px;" width="300" >
    </div>
Reveal cancer-associated immune cell states using unsupervised machine learning and statistical modeling. Derived from [Bachireddy et al. (2021) *Cell Reports*](https://doi.org/10.1016/j.celrep.2021.109992).

This repository has been newly refactored into a production-ready Python package, completely translating the original R-based methodology.

> **Note:** The original R scripts have been moved to the `legacy/` directory and are no longer maintained.

## Features

- **Common Factor Analysis (CFA)**: Maximum Likelihood Estimation (MLE) and varimax rotation to derive coherent gene expression modules.
- **Permutation Testing**: Evaluate empirical significance of factor correlations with specific gene sets.
- **Weighted T-Tests**: Handle imbalanced groups (e.g., Responder vs Non-Responder) with patient-adjusted weights and bootstrap FDR correction.
- **Meta-Clustering**: Pairwise Bhattacharyya distance metrics to merge and evaluate related cell clusters.

## Installation

Ensure you have Python 3.9+ installed.

```bash
git clone https://github.com/ozangocmen/cellstatemapping.git
cd cellstatemapping
pip install .
```

To install development dependencies (testing, formatting):
```bash
pip install .[dev]
```

## Quickstart (CLI)

The package provides a powerful Command Line Interface for direct execution without writing code.

**Run the full end-to-end pipeline:**
```bash
cellstatemapping run \
    --input data/patient1_counts.csv \
    --input data/patient2_counts.csv \
    --gene-set data/GS.proinf.txt \
    --metadata data/cell_metadata.csv \
    --groups Responder NonResponder \
    --output results_dir/
```

**Run only Common Factor Analysis:**
```bash
cellstatemapping cfa --gem expression.csv --gene-set GS.txt --n-factors 3
```

**Run only Weighted T-Tests:**
```bash
cellstatemapping test --metadata meta.csv --groups R NR --n-bootstrap 3000
```

## Python API Usage

You can also orchestrate the analysis using the robust Python API.

```python
from cellstatemapping import CellStateMappingPipeline

pipeline = CellStateMappingPipeline(output_dir="results/", n_factors=3)

# 1. Load data
pipeline.load_data(
    expression_files=["patient1.csv", "patient2.csv"],
    gene_set_file="GS.txt",
    metadata_file="metadata.csv"
)

# 2. Preprocess (log-normalize, select variable features)
pipeline.preprocess(n_variable_features=2000)

# 3. Factor Analysis
pipeline.run_factor_analysis()

# 4. Statistical significance tests
pipeline.run_permutation_tests(n_permutations=500)
pipeline.run_weighted_tests(group_names=("Responder", "NonResponder"))

# 5. Meta-clustering
pipeline.compute_metaclusters()

# Save everything
pipeline.save_results()
```

## Core Modules

If you'd like to use individual statistical components independently:

- `cellstatemapping.preprocessing`: Log-normalization (matches Seurat) and sampling
- `cellstatemapping.cfa`: `run_cfa`, `correlate_factors_with_geneset`
- `cellstatemapping.stats`: `weighted_t_test`, `compute_patient_weights`
- `cellstatemapping.distance`: `bhattacharyya_distance`, `compute_distance_matrix`

## Citation

```bibtex
@article{bachireddy2021mapping,
  title={Mapping the evolution of T cell states during response and resistance to adoptive cellular therapy},
  author={Bachireddy, P. and others},
  journal={Cell Reports},
  volume={37},
  number={6},
  pages={109992},
  year={2021},
  publisher={Elsevier}
}
```
