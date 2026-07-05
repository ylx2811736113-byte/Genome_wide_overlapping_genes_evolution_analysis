"""
Step 9. Build an ortholog matrix from species-specific ortholog assignment files.

Input:
    1. Species phylogeny in Newick format.
    2. Species directories containing final_ortholog_annotations.csv files.

Output:
    Ortholog matrix in CSV format, where each column
    represents a species and each row corresponds to
    ortholog assignment results.
"""

import os
import re
import pandas as pd


def extract_species_names(tree_file):
    """
    Extract species names from a Newick tree.
    """

    with open(tree_file, 'r') as file:
        tree_data = file.read()

    species_names = re.findall(
        r'([A-Za-z0-9_.]+)(?=[,:\)])',
        tree_data
    )

    return species_names


def build_ortholog_matrix(
        tree_file,
        root_folder,
        output_file):
    """
    Construct an ortholog matrix using species listed in the phylogenetic tree.
    """

    species_names = extract_species_names(tree_file)

    ortholog_matrix = pd.DataFrame()

    for species in species_names:

        species_folder = os.path.join(
            root_folder,
            species
        )

        ortholog_file = os.path.join(
            species_folder,
            'final_ortholog_annotations.csv'
        )

        if not os.path.exists(ortholog_file):

            print(
                f"Warning: ortholog file not found "
                f"for {species}"
            )
            continue

        df = pd.read_csv(ortholog_file)

        if 'Match' not in df.columns:

            print(
                f"Warning: 'Match' column not found "
                f"in {ortholog_file}"
            )
            continue

        # Add species-specific ortholog assignments
        species_data = df[['Match']].rename(
            columns={'Match': species}
        )

        ortholog_matrix = pd.concat(
            [ortholog_matrix, species_data],
            axis=1
        )

    ortholog_matrix.to_csv(
        output_file,
        index=False
    )

    print(
        f"Ortholog matrix saved to "
        f"{output_file}"
    )


if __name__ == "__main__":

    tree_file = (
        'D:/Work/overlapping_gene/revision/plant_tree.txt'
    )

    root_folder = (
        'D:/Work/overlapping_gene/revision/results/plant'
    )

    output_file = (
        'D:/Work/overlapping_gene/revision/results/plant_orthologues/result.csv'
    )

    build_ortholog_matrix(
        tree_file,
        root_folder,
        output_file
    )
