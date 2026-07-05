"""
Step 7. Finalize overlap annotations using adjacent-gene information.

Input:
    1. ortholog_overlap_refined.csv
    2. protein_with_position_update*.csv

Output:
    ortholog_overlap_final.csv containing the final
    overlap-state annotation for each
    orthologous gene pair.
"""

import os
import pandas as pd


# Root directory containing species-specific results
parent_folder = "D:/Work/overlapping_gene/results/animal_new"

processed_files = 0
skipped_files = 0


def check_positions(protein_string, position_dict):
    """
    Determine whether any protein in a semicolon-separated list is annotated as adjacent.
    """

    if pd.isna(protein_string):
        return "not adjacent"

    protein_string = str(protein_string).strip()

    if protein_string == "" or protein_string == "no":
        return "not adjacent"

    proteins = [
        p.strip()
        for p in protein_string.split(";")
        if p.strip() != "" and p.strip() != "no"
    ]

    positions = [
        position_dict.get(
            protein,
            "not adjacent"
        )
        for protein in proteins
    ]

    if "adjacent" in positions:
        return "adjacent"

    return "not adjacent"


# Process each species
for species in os.listdir(parent_folder):

    species_dir = os.path.join(
        parent_folder,
        species
    )

    if not os.path.isdir(species_dir):
        continue

    orthologs_file = os.path.join(
        species_dir,
        "ortholog_overlap_refined.csv"
    )

    position_file = next(
        (
            os.path.join(species_dir, f)
            for f in os.listdir(species_dir)
            if f.startswith(
                "protein_with_position_update"
            )
        ),
        None
    )

    output_file = os.path.join(
        species_dir,
        "ortholog_overlap_final.csv"
    )

    if not os.path.exists(orthologs_file):

        print(
            f"Skipping {species}: "
            f"ortholog_overlap_refined.csv not found."
        )

        skipped_files += 1
        continue

    if position_file is None:

        print(
            f"Skipping {species}: "
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
            f"Skipping {species}: "
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
            f"Skipping {species}: "
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

    ortholog_df["Protein1"] = (
        ortholog_df["Protein1"]
        .astype(str)
        .str.strip()
    )

    ortholog_df["Protein2"] = (
        ortholog_df["Protein2"]
        .astype(str)
        .str.strip()
    )

    # Refine pairs that are currently classified as Match = 0
    for index, row in ortholog_df.iterrows():

        if row["Match"] != 0:
            continue

        protein1 = row["Protein1"]
        protein2 = row["Protein2"]

        if protein1 == "no" or protein2 == "no":

            ortholog_df.at[
                index,
                "Match"
            ] = 0

            continue

        protein1_result = check_positions(
            protein1,
            protein_position_dict
        )

        protein2_result = check_positions(
            protein2,
            protein_position_dict
        )

        if (
            protein1_result == "adjacent"
            and
            protein2_result == "adjacent"
        ):

            ortholog_df.at[
                index,
                "Match"
            ] = 22

        elif (
            protein1_result == "adjacent"
            and
            protein2_result == "not adjacent"
        ):

            ortholog_df.at[
                index,
                "Match"
            ] = 21

        elif (
            protein1_result == "not adjacent"
            and
            protein2_result == "adjacent"
        ):

            ortholog_df.at[
                index,
                "Match"
            ] = 12

        else:

            ortholog_df.at[
                index,
                "Match"
            ] = 11

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
