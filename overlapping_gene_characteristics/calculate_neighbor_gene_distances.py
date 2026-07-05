import os
import pandas as pd

# Define the root directory
root_dir = "D:/Work/overlapping_gene/results/prokaryote_new"  # Replace with your actual root directory

# Traverse all subdirectories
for subdir, _, files in os.walk(root_dir):

    for file in files:

        # Process gene annotation files with neighboring gene information
        if file.endswith("extracted_genes_with_distances.csv"):

            file_path = os.path.join(subdir, file)

            print(f"Processing: {file_path}")

            try:
                # Read the annotation table
                df = pd.read_csv(file_path)

                # Check whether the required columns are present
                required_cols = {"chr_id", "Start", "End", "Strand"}

                if not required_cols.issubset(df.columns):
                    print(f"Missing required columns in {file_path}")
                    continue

                # Initialize new columns
                df["Gene_pre_len"] = pd.NA
                df["Gene_next_len"] = pd.NA
                df["Gene_pre_Strand"] = pd.NA
                df["Gene_next_Strand"] = pd.NA

                # Process each chromosome separately
                for chr_id, group in df.groupby("chr_id"):

                    group = group.sort_values("Start").copy()

                    for i in range(len(group)):

                        curr_idx = group.index[i]

                        current = group.iloc[i]

                        # Previous neighboring gene
                        if i > 0:

                            previous = group.iloc[i - 1]

                            df.loc[curr_idx, "Gene_pre_len"] = (
                                current["Start"] - previous["End"] - 1
                            )

                            df.loc[curr_idx, "Gene_pre_Strand"] = (
                                previous["Strand"] + current["Strand"]
                            )

                        # Next neighboring gene
                        if i < len(group) - 1:

                            following = group.iloc[i + 1]

                            df.loc[curr_idx, "Gene_next_len"] = (
                                following["Start"] - current["End"] - 1
                            )

                            df.loc[curr_idx, "Gene_next_Strand"] = (
                                current["Strand"] + following["Strand"]
                            )

                # Save the updated annotation table
                output_file = file_path.replace(
                    ".csv",
                    "_with_neighbor.csv"
                )

                df.to_csv(
                    output_file,
                    index=False,
                    encoding="utf-8"
                )

                print(f"Saved: {output_file}")

            except Exception as e:
                print(f"Error processing {file_path}: {e}")

print("Batch processing completed.")
