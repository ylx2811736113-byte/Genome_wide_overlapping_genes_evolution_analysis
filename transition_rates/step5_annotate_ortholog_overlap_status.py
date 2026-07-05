"""
Step 5. Annotate orthologous gene pairs according to overlap status.

Input:
    1. ortholog_pairs.csv
    2. overlapping_genes.csv

Output:
    ortholog_overlap_status.csv containing an additional
    'Match' column describing whether one or
    both orthologs overlap with genes in the
    target genome.
"""

import os
import pandas as pd


# Root directory containing species-specific results
main_folder = 'D:/Work/overlapping_gene/results/animal_new'

processed_files = 0
skipped_files = 0


def has_match(protein_string, protein_set):
    """
    Determine whether any protein ID in a semicolon-separated string is present in the overlapping protein set.
    """

    if pd.isna(protein_string):
        return False

    protein_string = str(protein_string).strip()

    if protein_string == '' or protein_string == 'no':
        return False

    proteins = protein_string.split(';')

    return any(
        protein in protein_set
        for protein in proteins
        if protein and protein != 'no'
    )


def process_species(species_dir):
    """
    Annotate ortholog pairs for one species.
    """

    global processed_files
    global skipped_files

    orthologs_file = os.path.join(
        species_dir,
        'ortholog_pairs.csv'
    )

    overlapping_file = os.path.join(
        species_dir,
        'overlapping_genes.csv'
    )

    if not os.path.exists(orthologs_file):

        print(
            f"Skipping {species_dir}: "
            f"ortholog_pairs.csv not found."
        )

        skipped_files += 1
        return

    if not os.path.exists(overlapping_file):

        print(
            f"Skipping {species_dir}: "
            f"overlapping_genes.csv not found."
        )

        skipped_files += 1
        return

    try:

        ortholog_df = pd.read_csv(
            orthologs_file,
            header=None,
            names=[
                'Gene1',
                'Protein1',
                'Gene2',
                'Protein2'
            ]
        )

        overlap_df = pd.read_csv(
            overlapping_file,
            header=None,
            names=['Protein']
        )

    except Exception as e:

        print(
            f"Skipping {species_dir}: "
            f"{e}"
        )

        skipped_files += 1
        return

    protein_set = set(
        overlap_df['Protein']
        .dropna()
        .astype(str)
    )

    ortholog_df['Match'] = 0

    for index, row in ortholog_df.iterrows():

        match1 = has_match(
            row['Protein1'],
            protein_set
        )

        match2 = has_match(
            row['Protein2'],
            protein_set
        )

        if match1 and match2:

            ortholog_df.at[
                index,
                'Match'
            ] = 33

        elif match1:

            ortholog_df.at[
                index,
                'Match'
            ] = 2

        elif match2:

            ortholog_df.at[
                index,
                'Match'
            ] = 4

    output_file = os.path.join(
        species_dir,
        'ortholog_overlap_status.csv'
    )

    ortholog_df.to_csv(
        output_file,
        index=False
    )

    processed_files += 1

    print(
        f"Saved: {output_file}"
    )


# Process all species
for species in os.listdir(main_folder):

    species_dir = os.path.join(
        main_folder,
        species
    )

    if os.path.isdir(species_dir):

        process_species(
            species_dir
        )

print("\nFinished.")

print(
    f"Files processed : {processed_files}"
)

print(
    f"Files skipped   : {skipped_files}"
)
