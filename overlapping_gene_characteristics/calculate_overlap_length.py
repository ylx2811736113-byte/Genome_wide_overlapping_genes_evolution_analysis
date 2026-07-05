import os
import pandas as pd

# Define the root directory
root_dir = "D:/Work/overlapping_gene/results/prokaryote_new"  # Replace with your actual root directory

# Traverse all subdirectories
for subdir, _, files in os.walk(root_dir):

    for file in files:

        # Process overlap annotation files
        if file.endswith("updated_overlap_genes.csv"):

            file_path = os.path.join(subdir, file)

            try:
                # Read the overlap annotation table
                df = pd.read_csv(file_path)

                # Calculate overlap length using the standard interval overlap formula
                df["Overlap_len"] = (
                    df[["Gene1_End", "Gene2_End"]].min(axis=1)
                    - df[["Gene1_Start", "Gene2_Start"]].max(axis=1)
                    + 1
                )

                # Save the updated table
                output_path = file_path.replace(
                    "updated_overlap_genes.csv",
                    "updated_overlap_genes_with_length.csv"
                )

                df.to_csv(output_path, index=False, encoding="utf-8")

            except Exception as e:
                print(f"Error processing {file_path}: {e}")

print("Batch processing completed.")
