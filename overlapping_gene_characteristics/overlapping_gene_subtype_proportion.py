import os
import pandas as pd

# Define the root directory and output file path
main_folder = "D:/Work/overlapping_gene/results/prokaryote_new/"  # Replace with your actual root directory
output_file = "D:/Work/overlapping_gene/results/prokaryote_new/output_summary.csv"

# Initialize a list to store the results
results = []

# Traverse all immediate subdirectories
for subfolder in os.listdir(main_folder):
    subfolder_path = os.path.join(main_folder, subfolder)

    # Process directories only
    if os.path.isdir(subfolder_path):

        # Traverse all files in the current subdirectory
        for file_name in os.listdir(subfolder_path):

            # Process files ending with "statistics.txt"
            if file_name.endswith("statistics.txt"):
                statistics_file = os.path.join(subfolder_path, file_name)

                # Skip if the file does not exist
                if os.path.isfile(statistics_file):

                    # Read the contents of the statistics file
                    with open(statistics_file, "r") as file:
                        lines = file.readlines()

                    # Extract overlap counts from the "Narrow_overlap_Type" section
                    overlap_data = {}
                    start_extracting = False

                    for line in lines:
                        line = line.strip()

                        if "Narrow_overlap_Type" in line:
                            start_extracting = True
                            continue

                        if start_extracting and line:
                            parts = line.split()
                            if len(parts) == 2:
                                category, count = parts[0], int(parts[1])
                                overlap_data[category] = count

                    # Ensure that all overlap categories are included
                    types = ["Co-oriented", "Convergent", "Divergent", "Nested"]
                    overlap_counts = {t: overlap_data.get(t, 0) for t in types}

                    # Calculate the total count and the proportion of each overlap category
                    total_count = sum(overlap_counts.values())
                    overlap_ratios = {
                        t: overlap_counts[t] / total_count if total_count > 0 else 0
                        for t in types
                    }

                    # Store the summary statistics for the current file
                    result = {
                        "Folder": subfolder,
                        "File": file_name,
                        "Co-oriented_Count": overlap_counts["Co-oriented"],
                        "Convergent_Count": overlap_counts["Convergent"],
                        "Divergent_Count": overlap_counts["Divergent"],
                        "Nested_Count": overlap_counts["Nested"],
                        "Co-oriented_Ratio": overlap_ratios["Co-oriented"],
                        "Convergent_Ratio": overlap_ratios["Convergent"],
                        "Divergent_Ratio": overlap_ratios["Divergent"],
                        "Nested_Ratio": overlap_ratios["Nested"],
                    }

                    # Append the result to the summary list
                    results.append(result)

# Save the summary table as a CSV file
results_df = pd.DataFrame(results)
results_df.to_csv(output_file, index=False)

print(f"Batch processing completed. Results have been saved to {output_file}")
