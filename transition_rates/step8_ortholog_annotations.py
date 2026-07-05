"""
Step 8. Finalize ortholog annotations.

This script scans all species directories, identifies
ortholog pairs containing missing ortholog assignments
(indicated by "no"), resets their Match value to 0,
and exports the finalized annotation table.

Input:
    ortholog_overlap_final.csv

Output:
    final_ortholog_annotations.csv
"""

import os
import pandas as pd


# Root directory containing species-specific result folders
main_folder = "D:/Work/overlapping_gene/results/animal_new"

processed_files = 0
skipped_files = 0


# Traverse all subdirectories
for root, _, files in os.walk(main_folder):

    if "ortholog_overlap_final.csv" not in files:
        continue

    input_file = os.path.join(
        root,
        "ortholog_overlap_final.csv"
    )

    output_file = os.path.join(
        root,
        "final_ortholog_annotations.csv"
    )

    try:

        df = pd.read_csv(input_file)

    except Exception as e:

        print(
            f"Skipping {input_file}: "
            f"{e}"
        )

        skipped_files += 1
        continue

    required_columns = [
        "Protein1",
        "Protein2",
        "Match"
    ]

    missing = [
        c for c in required_columns
        if c not in df.columns
    ]

    if missing:

        print(
            f"Skipping {input_file}: "
            f"missing columns "
            f"{', '.join(missing)}"
        )

        skipped_files += 1
        continue

    # Standardize missing-value representation
    df["Protein1"] = (
        df["Protein1"]
        .astype(str)
        .str.strip()
    )

    df["Protein2"] = (
        df["Protein2"]
        .astype(str)
        .str.strip()
    )

    # Reset Match for pairs lacking at least one ortholog
    missing_mask = (
        (df["Protein1"] == "no")
        |
        (df["Protein2"] == "no")
    )

    df.loc[
        missing_mask,
        "Match"
    ] = 0

    df.to_csv(
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
