"""
Step 6. Refine overlap annotations using adjacent-gene information.

Input:
    1. ortholog_overlap_status.csv
    2. protein_with_position_update*.csv

Output:
    ortholog_overlap_refined.csv with refined overlap
    classification based on neighboring
    gene information.
"""

import os
import pandas as pd


# Root directory containing species-specific results
base_folder = 'D:/Work/overlapping_gene/results/animal_new'

processed_files = 0
skipped_files = 0


def check_positions(protein_string, protein_position_dict):
    """
    Determine whether any protein in a semicolon-separated list is annotated as adjacent.
    """

    if pd.isna(protein_string):
        return "not adjacent"

    protein_string = str(protein_string).strip()

    if protein_string == "" or protein_string == "no":
        return "not adjacent"

    positions = []

    for protein in protein_string.split(";"):

        protein = protein.strip()

        if protein == "" or protein == "no":
            continue

        positions.append(
            protein_position_dict.get(
                protein,
                "not adjacent"
            )
        )

    if "adjacent" in positions:
        return "adjacent"

    return "not adjacent"


# Process each species
for folder_name in os.listdir(base_folder):

    folder_path = os.path.join(
        base_folder,
        folder_name
    )

    if not os.path.isdir(folder_path):
        continue

    orthologs_file = os.path.join(
        folder_path,
        "ortholog_overlap_status.csv"
    )

    position_file = next(
        (
            os.path.join(folder_path, f)
            for f in os.listdir(folder_path)
            if f.startswith(
                "protein_with_position_update"
            )
        ),
        None
    )

    if not os.path.exists(orthologs_file):

        print(
            f"Skipping {folder_name}: "
            f"ortholog_overlap_status.csv not found."
        )

        skipped_files += 1
        continue

    if position_file is None:

        print(
            f"Skipping {folder_name}: "
            f"position annotation file not found."
        )

        skipped_files += 1
        continue

    try:

        ortholog_df = pd.read_csv(
            orthologs_file
        )

        position_df = pd.read_csv(
            position_file
        )

    except Exception as e:

        print(
            f"Skipping {folder_name}: "
            f"{e}"
        )

        skipped_files += 1
        continue

    required_columns = [
        "Protein2",
        "Position"
    ]

    missing = [
        c for c in required_columns
        if c not in position_df.columns
    ]

    if missing:

        print(
            f"Skipping {folder_name}: "
            f"missing columns "
            f"{', '.join(missing)}"
        )

        skipped_files += 1
        continue

    protein_position_dict = (
        position_df
        .set_index("Protein2")["Position"]
        .to_dict()
    )

    for index, row in ortholog_df.iterrows():

        match_value = row["Match"]

        if match_value == 2:

            position = check_positions(
                row["Protein2"],
                protein_position_dict
            )

            ortholog_df.at[
                index,
                "Match"
            ] = 32 if position == "adjacent" else 31

        elif match_value == 4:

            position = check_positions(
                row["Protein1"],
                protein_position_dict
            )

            ortholog_df.at[
                index,
                "Match"
            ] = 23 if position == "adjacent" else 13

    output_file = os.path.join(
        folder_path,
        "ortholog_overlap_refined.csv"
    )

    ortholog_df.to_csv(
        output_file,
        index=False
    )

    processed_files += 1

    print(
        f"Saved: {output_file}"
    )

print("\nFinished.")

print(
    f"Files processed : {processed_files}"
)

print(
    f"Files skipped   : {skipped_files}"
)
