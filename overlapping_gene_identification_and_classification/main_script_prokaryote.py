"""
Pipeline for genome-wide identification and classification of overlapping protein-coding genes.

Species group:
    prokaryotes

Workflow:
    1. Extract CDS annotations.
    2. Generate CDS coordinate tables.
    3. Detect overlapping genes.
    4. Classify overlap patterns.

Input:
    GFF3 genome annotations.

Outputs:
    CDS coordinate table.
    Overlapping gene pairs.
    Overlap classification statistics.
"""

import os
from extract_cds_prokaryote import process_files as extract_cds
from identify_overlapping_genes import process_files as process_overlaps


def main(input_dir, output_dir):
    """Run the overlapping-gene identification workflow for all prokaryotic genome annotations."""

    # Scan genome directories
    folders = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]


    for folder in folders:
        folder_path = os.path.join(input_dir, folder)
        result_folder = os.path.join(output_dir, folder)

        # Create output directory
        if not os.path.exists(result_folder):
            os.makedirs(result_folder)

        # Locate annotation files
        gff_files = [f for f in os.listdir(folder_path) if f.endswith('.gff3') or f.endswith('.gff')]


        for gff_file in gff_files:

            gff_file_path = os.path.join(folder_path, gff_file)

            longest_cds_result_file = os.path.join(result_folder, f'{gff_file}_longest_cds.csv')

            print(f'Processing {gff_file_path}...')

            # Step 1. Extract CDS annotations
            extract_cds(gff_file_path, result_folder)

            # Step 2. Identify overlapping genes
            process_overlaps(longest_cds_result_file, result_folder)

            print(f'Finished processing {gff_file_path}\n')


if __name__ == "__main__":

    input_dir = '/public2/home/liulm/plant/plant_1'
    output_dir = '/public2/home/liulm/plant/result_1'


    main(input_dir, output_dir)
