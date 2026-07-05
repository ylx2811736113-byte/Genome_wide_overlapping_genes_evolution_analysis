import os
import pandas as pd

# Define the root directory and output file path
base_dir = "D:/Work/overlapping_gene/results/animal_new"  # Replace with your actual root directory
output_file = os.path.join(base_dir, "remove_intron_length_summary.csv")

# Initialize a list to store the results
results = []

# Traverse all immediate subdirectories
for folder_name in os.listdir(base_dir):

    folder_path = os.path.join(base_dir, folder_name)

    # Process directories only
    if not os.path.isdir(folder_path):
        continue

    input_file = os.path.join(
        folder_path,
        "genes_remove_intron_length.csv"
    )

    # Skip if the input file does not exist
    if not os.path.exists(input_file):
        print(f"File not found: {input_file}")
        continue

    try:
        # Read the input table
        df = pd.read_csv(input_file)

        # Check whether the required columns are present
        required_cols = {
            "og_remove_intron_length",
            "nonog_remove_intron_length",
        }

        if not required_cols.issubset(df.columns):
            print(f"Missing required columns in {input_file}")
            continue

        # Calculate the average intron-removed length of overlapping genes
        og_mean = (
            df.loc[
                df["og_remove_intron_length"] > 0,
                "og_remove_intron_length",
            ]
            .mean()
        )

        # Calculate the average intron-removed length of non-overlapping genes
        nonog_mean = (
            df.loc[
                df["nonog_remove_intron_length"] > 0,
                "nonog_remove_intron_length",
            ]
            .mean()
        )

        # Store the summary statistics
        results.append(
            {
                "Species": folder_name,
                "Average_OG_RemoveIntron_Length": og_mean,
                "Average_NonOG_RemoveIntron_Length": nonog_mean,
            }
        )

    except Exception as e:
        print(f"Error processing {input_file}: {e}")

# Save the summary table as a CSV file
result_df = pd.DataFrame(results)

result_df.to_csv(
    output_file,
    index=False,
    encoding="utf-8"
)

print(f"Results have been saved to {output_file}")
