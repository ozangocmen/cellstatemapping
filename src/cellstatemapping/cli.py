"""
Command-line interface for cellstatemapping.

Provides CLI commands for running the cell state mapping pipeline
from the terminal.

Usage
-----
    cellstatemapping run --input data/ --gene-set GS.txt --output results/
    cellstatemapping cfa --gem gem.csv --gene-set GS.txt --n-factors 3
    cellstatemapping test --metadata metadata.csv --groups R NR
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import click

from cellstatemapping.pipeline import CellStateMappingPipeline

logger = logging.getLogger(__name__)


def _setup_logging(verbose: bool) -> None:
    """Configure logging based on verbosity."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        stream=sys.stderr,
    )


@click.group()
@click.version_option(version="1.0.0", prog_name="cellstatemapping")
def main():
    """Cell State Mapping — Reveal cancer-associated immune cell states.

    A Python toolkit for unsupervised machine learning and statistical
    modeling of single-cell immune cell states.
    """
    pass


@main.command()
@click.option(
    "--input", "-i",
    "input_files",
    required=True,
    multiple=True,
    type=click.Path(exists=True),
    help="Input expression matrix CSV file(s).",
)
@click.option(
    "--gene-set", "-g",
    "gene_set_file",
    required=True,
    type=click.Path(exists=True),
    help="Gene set file (one gene per line).",
)
@click.option(
    "--metadata", "-m",
    "metadata_file",
    required=False,
    type=click.Path(exists=True),
    help="Metadata CSV file with cluster, patient, group columns.",
)
@click.option(
    "--output", "-o",
    "output_dir",
    default="results",
    type=click.Path(),
    help="Output directory for results.",
)
@click.option("--n-factors", default=3, type=int, help="Number of CFA factors.")
@click.option("--n-permutations", default=500, type=int, help="Permutation test iterations.")
@click.option("--n-bootstrap", default=3000, type=int, help="Bootstrap iterations for p-value correction.")
@click.option("--groups", nargs=2, default=("R", "NR"), help="Group labels (responder non-responder).")
@click.option("--seed", default=42, type=int, help="Random seed.")
@click.option("--verbose", "-v", is_flag=True, help="Enable debug logging.")
def run(
    input_files,
    gene_set_file,
    metadata_file,
    output_dir,
    n_factors,
    n_permutations,
    n_bootstrap,
    groups,
    seed,
    verbose,
):
    """Run the full cell state mapping pipeline."""
    _setup_logging(verbose)

    click.echo(f"Cell State Mapping Pipeline v1.0.0")
    click.echo(f"Input files: {input_files}")
    click.echo(f"Gene set: {gene_set_file}")
    click.echo(f"Output: {output_dir}")

    pipeline = CellStateMappingPipeline(
        output_dir=output_dir,
        n_factors=n_factors,
        n_permutations=n_permutations,
        n_bootstrap=n_bootstrap,
        seed=seed,
    )

    # Load data
    pipeline.load_data(
        expression_files=list(input_files),
        gene_set_file=gene_set_file,
        metadata_file=metadata_file,
    )

    # Preprocess
    pipeline.preprocess()

    # Factor analysis
    click.echo("Running factor analysis...")
    pipeline.run_factor_analysis()

    # Permutation tests
    click.echo("Running permutation tests...")
    pipeline.run_permutation_tests()

    # Weighted tests (if metadata provided)
    if metadata_file:
        click.echo("Running weighted cluster comparisons...")
        pipeline.run_weighted_tests(group_names=tuple(groups))

    # Meta-clustering
    click.echo("Computing meta-clusters...")
    pipeline.compute_metaclusters()

    # Save
    pipeline.save_results()
    click.echo(f"Results saved to {output_dir}/")


@main.command()
@click.option("--gem", required=True, type=click.Path(exists=True), help="Expression matrix CSV.")
@click.option("--gene-set", required=True, type=click.Path(exists=True), help="Gene set file.")
@click.option("--n-factors", default=3, type=int, help="Number of factors.")
@click.option("--output", "-o", default="results", type=click.Path(), help="Output directory.")
@click.option("--verbose", "-v", is_flag=True, help="Enable debug logging.")
def cfa(gem, gene_set, n_factors, output, verbose):
    """Run only Common Factor Analysis."""
    _setup_logging(verbose)

    from cellstatemapping.preprocessing import load_expression_matrix, load_gene_set
    from cellstatemapping.cfa import correlate_factors_with_geneset

    click.echo(f"Running CFA with {n_factors} factors...")

    gem_df = load_expression_matrix(gem, transpose=True).T  # cells x genes
    gs = load_gene_set(gene_set)

    loadings, correlations = correlate_factors_with_geneset(
        gem_df, gs, n_factors=n_factors
    )

    output_dir = Path(output)
    output_dir.mkdir(parents=True, exist_ok=True)
    loadings.to_csv(output_dir / "factor_loadings.csv")

    click.echo("Factor–Gene Set Correlations:")
    for factor, (r, p) in correlations.items():
        click.echo(f"  {factor}: r={r:.4f}, p={p:.6f}")

    click.echo(f"Loadings saved to {output_dir}/factor_loadings.csv")


@main.command()
@click.option("--metadata", required=True, type=click.Path(exists=True), help="Metadata CSV.")
@click.option("--groups", nargs=2, default=("R", "NR"), help="Group labels.")
@click.option("--n-bootstrap", default=3000, type=int, help="Bootstrap iterations.")
@click.option("--output", "-o", default="results", type=click.Path(), help="Output directory.")
@click.option("--verbose", "-v", is_flag=True, help="Enable debug logging.")
def test(metadata, groups, n_bootstrap, output, verbose):
    """Run only weighted t-test comparison."""
    _setup_logging(verbose)
    import pandas as pd
    from cellstatemapping.stats import weighted_cluster_comparison

    click.echo(f"Running weighted t-tests: {groups[0]} vs {groups[1]}...")

    meta_df = pd.read_csv(metadata, index_col=0)

    results = weighted_cluster_comparison(
        metadata=meta_df,
        group_names=tuple(groups),
        n_bootstrap=n_bootstrap,
    )

    output_dir = Path(output)
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "weighted_test_results.csv", index=False)

    click.echo("Results:")
    click.echo(results.to_string(index=False))
    click.echo(f"\nSaved to {output_dir}/weighted_test_results.csv")


if __name__ == "__main__":
    main()
