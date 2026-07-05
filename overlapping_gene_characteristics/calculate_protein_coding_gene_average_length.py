import os
import pandas as pd

# Define the root directory and output file path
root_dir = "D:/Work/overlapping_gene/results/prokaryote_new"  # Replace with your actual root directory
output_file = "D:/Work/overlapping_gene/results/prokaryote_new/protein_coding_gene_average_length.csv"

# Initialize a dictionary to store the results
results = {}

# Traverse all subdirectories
for subdir, _, files in os.walk(root_dir):

    for file in files:

        # Process gene annotation files
        if file.endswith("extracted_genes.csv"):

            file_path = os.path.join(subdir, file)

            try:
                # Read the gene annotation table
                df = pd.read_csv(file_path, sep="\t")

                # Calculate gene lengths if the required columns are present
                if {"Start", "End"}.issubset(df.columns):

                    df["Length"] = (
                        df["End"] - df["Start"]
                    ).abs() + 1

                    avg_length = df["Length"].mean()

                else:
                    print(f"Missing required columns in {file_path}")
                    avg_length = float("nan")

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
    columns=["Average_ProteinCoding_Gene_Length"],
)

result_df.index.name = "subfolder"

result_df.to_csv(
    output_file,
    encoding="utf-8"
)

print(f"Results have been saved to {output_file}")
