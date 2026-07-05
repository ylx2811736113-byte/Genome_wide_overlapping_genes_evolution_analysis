"""
Step 10. Generate species trait files from the ortholog matrix.

Input:
    Ortholog matrix generated in Step 9.

Output:
    One trait file per ortholog group, where rows
    correspond to species and columns represent
    binary trait states for downstream phylogenetic
    analyses.
"""

import os
import pandas as pd


def generate_trait_files(
        ortholog_matrix_file,
        output_dir):
    """
    Convert an ortholog matrix into individual trait files for phylogenetic analyses.
    """

    os.makedirs(
        output_dir,
        exist_ok=True
    )

    trait_matrix = pd.read_csv(
        ortholog_matrix_file,
        sep=","
    )

    if "Species" not in trait_matrix.columns:
        raise ValueError(
            "Input matrix must contain a 'Species' column."
        )

    # Generate one trait file per ortholog group
    for ortholog_group in trait_matrix.columns[1:]:

        trait_data = trait_matrix[
            ["Species", ortholog_group]
        ]

        output_file = os.path.join(
            output_dir,
            f"{ortholog_group}_traits.txt"
        )

        trait_data.to_csv(
            output_file,
            sep="\t",
            index=False,
            header=False
        )

        print(
            f"Trait file generated: "
            f"{output_file}"
        )


if __name__ == "__main__":

    ortholog_matrix_file = (
        "D:/Work/overlapping_gene/revision/results/plant_orthologues/result.csv"
    )

    output_dir = (
        "D:/Work/overlapping_gene/revision/results/traits_files_plant_020426"
    )

    generate_trait_files(
        ortholog_matrix_file,
        output_dir
    )
