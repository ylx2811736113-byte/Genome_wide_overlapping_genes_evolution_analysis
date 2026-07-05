"""
Calculate the tissue specificity index (Tau) for each gene in paired gene expression data.

Input:
    A CSV file containing paired gene expression values where:
        - Expression columns for Gene 1 end with "_g1".
        - Expression columns for Gene 2 end with "_g2".

Output:
    The input table with two additional columns:
        - Tau1: Tau index for Gene 1
        - Tau2: Tau index for Gene 2

Reference:
    Yanai et al. (2005)
    Genome-wide midrange transcription profiles reveal expression level relationships in human tissue specification.
"""

from pathlib import Path
import numpy as np
import pandas as pd


def tau_index(expression_values):
    """
    Calculate the tissue specificity index (Tau).

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
    """Calculate Tau indices for paired genes."""

    input_file = Path(
        "D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/random_control_update.csv"
    )

    output_file = Path(
        "D:/Work/overlapping_gene/Rfiles/gene_expression/human_061526_1/tau_random_control.csv"
    )

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found:\n{input_file}")

    # Read input data
    df = pd.read_csv(input_file)

    if df.empty:
        raise ValueError("The input file is empty.")

    # Identify expression columns
    g1_cols = [col for col in df.columns if col.endswith("_g1")]
    g2_cols = [col for col in df.columns if col.endswith("_g2")]

    if not g1_cols:
        raise ValueError("No Gene 1 expression columns ending with '_g1' were found.")

    if not g2_cols:
        raise ValueError("No Gene 2 expression columns ending with '_g2' were found.")

    if len(g1_cols) != len(g2_cols):
        raise ValueError(
            "The numbers of Gene 1 and Gene 2 expression columns do not match."
        )

    # Calculate Tau indices
    df["Tau1"] = df[g1_cols].apply(tau_index, axis=1)
    df["Tau2"] = df[g2_cols].apply(tau_index, axis=1)

    # Save results
    df.to_csv(output_file, index=False)

    print("Tau index calculation completed successfully.")
    print(f"Results saved to:\n{output_file}")


if __name__ == "__main__":
    main()
