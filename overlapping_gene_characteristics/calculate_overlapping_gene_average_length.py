import os
import pandas as pd

# Define the root directory and output file path
root_dir = "D:/Work/overlapping_gene/results/prokaryote_new"  # Replace with your actual root directory
output_file = "D:/Work/overlapping_gene/results/prokaryote_new/average_overlap_gene_lengths.csv"

# Initialize a dictionary to store the results
results = {}

# Traverse all subdirectories
for subdir, _, files in os.walk(root_dir):

    for file in files:
        if file.endswith("updated_overlap_genes.csv"):

            file_path = os.path.join(subdir, file)

            try:
                # Read overlapping gene pairs
                df = pd.read_csv(file_path, sep=",")

                # Extract genomic coordinates for Gene1
                gene1_df = df[["Gene1_ID", "Gene1_Start", "Gene1_End"]].copy()
                gene1_df.columns = ["Gene_ID", "Start", "End"]

                # Extract genomic coordinates for Gene2
                gene2_df = df[["Gene2_ID", "Gene2_Start", "Gene2_End"]].copy()
                gene2_df.columns = ["Gene_ID", "Start", "End"]

                # Combine the two gene sets
                combined_df = pd.concat([gene1_df, gene2_df], ignore_index=True)

                # Remove duplicated genes while preserving a copy of the data
                unique_genes = combined_df.drop_duplicates(
                    subset=["Gene_ID"]
                ).copy()

                # Calculate gene lengths and their average
                unique_genes["Length"] = (
                    abs(unique_genes["End"] - unique_genes["Start"]) + 1
                )
                avg_length = unique_genes["Length"].mean()

            except Exception as e:
                print(f"Error processing {file_path}: {e}")
                avg_length = float("nan")

            # Store the average length using the subdirectory name as the key
            folder_name = os.path.basename(subdir)
            results[folder_name] = avg_length

# Save the summary table as a CSV file
result_df = pd.DataFrame.from_dict(
    results,
    orient="index",
    columns=["Average_Overlapping_Gene_Length"],
)
result_df.index.name = "subfolder"
result_df.to_csv(output_file)

print(f"Results have been saved to {output_file}")
