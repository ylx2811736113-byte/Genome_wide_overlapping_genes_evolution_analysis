"""
Step 1. Calculate intergenic distances between adjacent genes.

Input:
    Multiple extracted_genes.csv files.

Output:
    Updated gene tables with an additional
    'distance_to_next' column representing the
    genomic distance to the next gene located on
    the same chromosome or scaffold.
"""

import os
import glob
import pandas as pd


# Root directory containing genome-specific results
root_dir = 'D:/Work/overlapping_gene/results/animal_new'

# Locate all extracted gene tables
file_pattern = os.path.join(
    root_dir,
    '**',
    '*extracted_genes.csv'
)

file_list = glob.glob(
    file_pattern,
    recursive=True
)

if len(file_list) == 0:
    raise FileNotFoundError(
        "No extracted_genes.csv files were found."
    )

processed_files = 0
skipped_files = 0

# Process each genome
for file_path in file_list:

    try:
        df = pd.read_csv(file_path)

    except Exception as e:
        print(
            f"Skipping {file_path}: "
            f"unable to read file ({e})"
        )
        skipped_files += 1
        continue

    # Check required columns
    required_columns = [
        'chr_id',
        'Start',
        'End'
    ]

    missing_columns = [
        col for col in required_columns
        if col not in df.columns
    ]

    if missing_columns:
        print(
            f"Skipping {file_path}: "
            f"missing columns "
            f"{', '.join(missing_columns)}."
        )
        skipped_files += 1
        continue

    if df.empty:
        print(
            f"Skipping {file_path}: "
            f"empty file."
        )
        skipped_files += 1
        continue

    # Sort genes by genomic coordinates
    df = (
        df.sort_values(
            by=['chr_id', 'Start']
        )
        .reset_index(drop=True)
    )

    # Coordinates of the next gene
    next_chr = df['chr_id'].shift(-1)
    next_start = df['Start'].shift(-1)

    # Calculate intergenic distances
    df['distance_to_next'] = (
        next_start - df['End']
    ).where(
        df['chr_id'] == next_chr
    )

    output_file = file_path.replace(
        'extracted_genes.csv',
        'extracted_genes_with_distances.csv'
    )

    df.to_csv(
        output_file,
        index=False
    )

    processed_files += 1

    print(
        f"Processed genome: "
        f"{os.path.basename(os.path.dirname(file_path))}"
    )

print("\nFinished.")

print(
    f"Processed files : {processed_files}"
)

print(
    f"Skipped files   : {skipped_files}"
)
