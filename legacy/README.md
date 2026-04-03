# Legacy Scripts (Deprecated)

> **⚠️ These scripts are deprecated and preserved for historical reference only.**

This directory contains the original R scripts and R Markdown files that were used during the
initial development of the `cellstatemapping` methodology, derived from
[Bachireddy et al. (2021) Cell Reports](https://doi.org/10.1016/j.celrep.2021.109992).

These scripts are **not maintained** and should **not be used** for new analyses.
Please use the Python `cellstatemapping` package instead:

```python
from cellstatemapping import CellStateMappingPipeline
```

## Directory Structure

- `Objective1/` — Non-PCA-based factor reproduction, annotation, and interpretation
  - `BIGQA-WP1-functions.R` — Core R functions for CFA and permutation testing
  - `CFA-Stats.R` — CFA statistics and exploratory analysis
  - `automatization.R` — Automated data loading for multiple patient samples
- `Objective2/` — Weighted one-sided t-tests and statistical modeling
  - `weighted_t_test.R` — Weighted t-test implementation with bootstrap correction
  - `wtd_t_test.Rmd` — Documentation and examples for weighted t-test
  - `statsmodel.ipynb` — Alternative Python implementation (now integrated into the package)
- `Objective3/` — Meta-clustering via Bhattacharyya distance
  - `Objective 3.Rmd` — Bhattacharyya distance computation and heatmap visualization
- `final_script.Rmd` — End-to-end pipeline combining all objectives

## Citation

If you use this software, please cite:

> Bachireddy, P., et al. (2021). Mapping the evolution of T cell states during response and
> resistance to adoptive cellular therapy. *Cell Reports*, 37(6), 109992.
