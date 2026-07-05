import os
import pandas as pd

# Define the root directory and output file path
root_dir = "D:/Work/overlapping_gene/results/prokaryote_new/"  # Replace with your actual root directory
output_file = "D:/Work/overlapping_gene/results/prokaryote_new/ogs_summary.csv"

# Initialize a list to store the results
results = []

# Traverse all subdirectories
for folder_name, subfolders, filenames in os.walk(root_dir):
    for filename in filenames:
        if filename.startswith("chr_") and filename.endswith(".csv"):
            file_path = os.path.join(folder_name, filename)

            try:
                # Read the CSV file and calculate the sum of the second column
                # (i.e., "the number of ogs")
                df = pd.read_csv(file_path)
                ogs_sum = df.iloc[:, 1].sum()

                # Store the subdirectory name and the corresponding total OGS count
                results.append([os.path.basename(folder_name), ogs_sum])

            except Exception as e:
                print(f"Error processing file {file_path}: {e}")

# Save the summary table as a CSV file
df_result = pd.DataFrame(results, columns=["Folder_Name", "OGS_Sum"])
df_result.to_csv(output_file, index=False, encoding="utf-8")

print(f"Results have been saved to {output_file}")
