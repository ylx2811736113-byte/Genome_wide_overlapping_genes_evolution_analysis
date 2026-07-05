"""
Calculate the tissue specificity index (Tau) for each gene based on its expression levels across multiple tissues.

Input:
    A CSV file where:
        - The first column contains gene identifiers.
        - All remaining columns contain gene expression values for different tissues.

Output:
    The input table with an additional column named "Tau".

Reference:
    Yanai et al. (2005)
    Genome-wide midrange transcription profiles reveal expression level relationships in human tissue specification.
"""

from pathlib import Path
import numpy as np
import pandas as pd


def tau_index(expression_values):
    """
    Calculate the tissue specificity index (Tau) for one gene.

    Parameters
    ----------
    expression_values : array-like
        Gene expression values across tissues.

    Returns
    -------
    float
        Tau index ranging from 0 to 1.
        - 0 indicates ubiquitous expression.
        - 1 indicates tissue-specific expression.
        Returns NaN if all expression values are zero or missing.
    """
    expr = pd.to_numeric(expression_values, errors="coerce").to_numpy(dtype=float)

    # Replace missing values with zero
    expr = np.nan_to_num(expr, nan=0.0)

    max_expr = np.max(expr)

    # Avoid division by zero
    if max_expr == 0:
        return np.nan

    # At least two tissues are required
    if len(expr) < 2:
        return np.nan

    normalized_expr = expr / max_expr

    tau = np.sum(1 - normalized_expr) / (len(normalized_expr) - 1)

    return tau


def main():
    """Calculate Tau index for all genes."""

    input_file = Path(
        "D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/all_genes_simple_update.csv"
    )

    output_file = Path(
        "D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/all_genes.csv"
    )

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found:\n{input_file}")

    # Read expression matrix
    df = pd.read_csv(input_file)

    if df.empty:
        raise ValueError("The input file is empty.")

    if df.shape[1] < 2:
        raise ValueError(
            "The input file must contain at least one gene identifier column and one expression column."
        )

    # Expression columns (all columns except the first one)
    expression_matrix = df.iloc[:, 1:]

    # Calculate Tau for each gene
    df["Tau"] = expression_matrix.apply(tau_index, axis=1)

    # Save results
    df.to_csv(output_file, index=False)

    print(f"Tau index calculation completed successfully.")
    print(f"Results saved to:\n{output_file}")


if __name__ == "__main__":
    main()
