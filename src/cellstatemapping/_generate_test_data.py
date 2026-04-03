"""Generate minimal synthetic test data for cellstatemapping tests.

Run this script to regenerate test data files in src/cellstatemapping/data/:
    py -m cellstatemapping._generate_test_data
"""

import numpy as np
import pandas as pd
from pathlib import Path


def generate_test_data():
    """Generate synthetic expression matrix, metadata, and gene set files."""
    data_dir = Path(__file__).parent / "data"
    data_dir.mkdir(exist_ok=True)

    rng = np.random.default_rng(42)
    n_cells = 200
    n_genes = 100
    n_patients = 6
    n_clusters = 4

    # Gene names
    gene_names = [f"Gene{i}" for i in range(1, n_genes + 1)]
    cell_names = [f"Cell{i}" for i in range(1, n_cells + 1)]

    # Expression matrix: cells x genes (Poisson-distributed counts)
    counts = rng.poisson(lam=5, size=(n_cells, n_genes))
    gem = pd.DataFrame(counts, index=cell_names, columns=gene_names)
    gem.to_csv(data_dir / "test_expression_matrix.csv")

    # Gene set (subset of genes)
    gene_set = gene_names[:20]
    with open(data_dir / "test_gene_set.txt", "w") as f:
        for g in gene_set:
            f.write(g + "\n")

    # Curated gene sets
    with open(data_dir / "test_curated_gene_sets.txt", "w") as f:
        f.write("GeneSet1\tTest gene set 1\t" + "\t".join(gene_set[:10]) + "\n")
        f.write("GeneSet2\tTest gene set 2\t" + "\t".join(gene_set[10:20]) + "\n")

    # Metadata
    patients = [f"Patient{rng.integers(1, n_patients + 1)}" for _ in range(n_cells)]
    groups = []
    for p in patients:
        pid = int(p.replace("Patient", ""))
        groups.append("R" if pid <= 3 else "NR")
    clusters = [f"Cluster{rng.integers(1, n_clusters + 1)}" for _ in range(n_cells)]

    metadata = pd.DataFrame({
        "patient": patients,
        "group": groups,
        "cluster": clusters,
    }, index=cell_names)
    metadata.to_csv(data_dir / "test_metadata.csv")

    print(f"Generated test data in {data_dir}:")
    print(f"  - test_expression_matrix.csv ({n_cells} cells x {n_genes} genes)")
    print(f"  - test_gene_set.txt ({len(gene_set)} genes)")
    print(f"  - test_curated_gene_sets.txt (2 gene sets)")
    print(f"  - test_metadata.csv ({n_cells} cells, {n_patients} patients)")


if __name__ == "__main__":
    generate_test_data()
