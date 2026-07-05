import os
import pandas as pd

# Define the root directory and output file path
root_dir = "D:/Work/overlapping_gene/results/prokaryote_new/"  # Replace with your actual root directory
output_file = "D:/Work/overlapping_gene/results/prokaryote_new/protein_coding_genes_count.csv"

# Initialize a list to store the results
results = []

# Traverse all subdirectories
for folder_name, subfolders, filenames in os.walk(root_dir):
    for filename in filenames:
        if filename.endswith("extracted_genes.csv"):
            file_path = os.path.join(folder_name, filename)

            try:
                # Count the number of data rows (excluding the header)
                total_lines = sum(1 for _ in open(file_path)) - 1

                # Store the subdirectory name and the corresponding gene count
                results.append([os.path.basename(folder_name), total_lines])

            except Exception as e:
                print(f"Error processing file {file_path}: {e}")

# Save the summary table as a CSV file
df = pd.DataFrame(results, columns=["Folder_Name", "Row_Count"])
df.to_csv(output_file, index=False, encoding="utf-8")

print(f"Results have been saved to {output_file}")
