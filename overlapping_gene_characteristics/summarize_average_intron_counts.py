import os
import pandas as pd

# Define the root directory and output file path
base_dir = "D:/Work/overlapping_gene/results/animal_new"  # Replace with your actual root directory
output_file = os.path.join(base_dir, "intron_counts_summary.csv")

# Initialize a list to store the results
results = []

# Traverse all immediate subdirectories
for folder_name in os.listdir(base_dir):

    folder_path = os.path.join(base_dir, folder_name)

    # Process directories only
    if not os.path.isdir(folder_path):
        continue

    # Define the input file paths
    og_file = os.path.join(folder_path, "og_matched_genes_intron_with_counts.csv")
    nonog_file = os.path.join(folder_path, "og_unmatched_genes_intron_with_counts.csv")

    # Skip if either input file is missing
    if not (os.path.exists(og_file) and os.path.exists(nonog_file)):
        print(f"Missing input files in {folder_name}")
        continue

    try:
        # Read the input tables
        og_df = pd.read_csv(og_file)
        nonog_df = pd.read_csv(nonog_file)

        # Check whether the required column is present
        required_col = "intron_counts"

        if (
            required_col not in og_df.columns
            or required_col not in nonog_df.columns
        ):
            print(f"Missing required column in {folder_name}")
            continue

        # Calculate the average intron count of overlapping genes
        og_mean = og_df["intron_counts"].mean()

        # Calculate the average intron count of non-overlapping genes
        nonog_mean = nonog_df["intron_counts"].mean()

        # Store the summary statistics
        results.append(
            {
                "Species": folder_name,
                "Average_OG_Intron_Count": og_mean,
                "Average_NonOG_Intron_Count": nonog_mean,
            }
        )

    except Exception as e:
        print(f"Error processing {folder_name}: {e}")

# Save the summary table as a CSV file
result_df = pd.DataFrame(results)

result_df.to_csv(
    output_file,
    index=False,
    encoding="utf-8",
)

print(f"Results have been saved to {output_file}")
