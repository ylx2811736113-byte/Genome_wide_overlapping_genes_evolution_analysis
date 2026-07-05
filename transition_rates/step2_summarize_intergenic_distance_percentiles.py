"""
Step 2. Summarize intergenic distance statistics across genomes.

Input:
    Multiple extracted_genes_with_distances.csv files.

Output:
    Summary table containing the 5th percentile,
    10th percentile, and mean positive intergenic
    distance for each genome.
"""

import os
import glob
import math
import pandas as pd


# Root directory containing genome-specific results
root_dir = 'D:/Work/overlapping_gene/results/animal_new'

# Store summary statistics
summary_data = []

# Locate all intergenic distance tables
file_pattern = os.path.join(
    root_dir,
    '**',
    '*extracted_genes_with_distances.csv'
)

file_list = glob.glob(
    file_pattern,
    recursive=True
)

if len(file_list) == 0:
    raise FileNotFoundError(
        "No extracted_genes_with_distances.csv files were found."
    )

# Process each genome
for file_path in file_list:

    folder_name = os.path.basename(
        os.path.dirname(file_path)
    )

    try:
        df = pd.read_csv(file_path)

    except Exception as e:
        print(
            f"Skipping {folder_name}: "
            f"unable to read file ({e})"
        )
        continue

    # Check required column
    if 'distance_to_next' not in df.columns:
        print(
            f"Skipping {folder_name}: "
            f"'distance_to_next' column not found."
        )
        continue

    # Retain positive intergenic distances
    positive_distances = df.loc[
        df['distance_to_next'] > 0,
        'distance_to_next'
    ]

    if positive_distances.empty:
        print(
            f"Skipping {folder_name}: "
            f"no positive intergenic distances."
        )
        continue

    sorted_distances = (
        positive_distances
        .sort_values()
        .reset_index(drop=True)
    )

    n = len(sorted_distances)

    # Determine percentile positions
    index_5 = max(
        0,
        min(math.ceil(n * 0.05) - 1, n - 1)
    )

    index_10 = max(
        0,
        min(math.ceil(n * 0.10) - 1, n - 1)
    )

    value_5 = sorted_distances.iloc[index_5]
    value_10 = sorted_distances.iloc[index_10]

    mean_distance = positive_distances.mean()

    summary_data.append({

        'Folder Name': folder_name,

        '5% Value': value_5,
        '5% Index': index_5 + 1,

        '10% Value': value_10,
        '10% Index': index_10 + 1,

        'Mean Value (>0)': mean_distance

    })

    print(
        f"Processed genome: {folder_name}"
    )

# Export summary table
if summary_data:

    summary_df = pd.DataFrame(summary_data)

    output_file = os.path.join(
        root_dir,
        'intergenic_distance_thresholds.csv'
    )

    summary_df.to_csv(
        output_file,
        index=False
    )

    print(
        f"\nSummary table saved to:\n{output_file}"
    )

else:

    print(
        "\nNo valid genomes were processed."
    )
