import os
import pandas as pd

# Define the root directory and output file path
root_dir = "D:/Work/overlapping_gene/results/prokaryote_new"  # Replace with your actual root directory
output_file = "D:/Work/overlapping_gene/results/prokaryote_new/non_overlapping_gene_length_avg.csv"

# Initialize a dictionary to store the results
results = {}

# Traverse all subdirectories
for subdir, _, files in os.walk(root_dir):

    # Locate the required input files
    extracted_file = None
    overlap_file = None

    for file in files:
        if file.endswith("extracted_genes.csv"):
            extracted_file = os.path.join(subdir, file)
        elif file.endswith("updated_overlap_genes.csv"):
            overlap_file = os.path.join(subdir, file)

    # Process the data only if both input files are available
    if extracted_file and overlap_file:

        try:
            # Read gene annotation data
            df_extracted = pd.read_csv(extracted_file, sep="\t")
            df_extracted = df_extracted[["GeneID", "Start", "End"]].copy()

            # Read overlapping gene pairs
            df_overlap = pd.read_csv(overlap_file, sep=",")
            overlap_gene_ids = set(df_overlap["Gene1_ID"]) | set(df_overlap["Gene2_ID"])

            # Retain only non-overlapping genes
            df_non_overlap = df_extracted[
                ~df_extracted["GeneID"].isin(overlap_gene_ids)
            ].copy()

            # Calculate gene lengths and their average
            df_non_overlap["Length"] = (
                abs(df_non_overlap["End"] - df_non_overlap["Start"]) + 1
            )
            avg_length = df_non_overlap["Length"].mean()

        except Exception as e:
            print(f"Error processing {subdir}: {e}")
            avg_length = float("nan")

        # Store the average length using the subdirectory name as the key
        folder_name = os.path.basename(subdir)
        results[folder_name] = avg_length

# Save the summary table as a CSV file
result_df = pd.DataFrame.from_dict(
    results,
    orient="index",
    columns=["Average_NonOverlapping_Gene_Length"],
)
result_df.index.name = "subfolder"
result_df.to_csv(output_file)

print(f"Results have been saved to {output_file}")
