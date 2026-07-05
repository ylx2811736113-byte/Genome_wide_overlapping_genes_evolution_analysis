"""
Pipeline for genome-wide identification and classification of overlapping protein-coding genes.

Species group:
    Animal

Workflow:
    1. Parse genome annotations.
    2. Select longest transcripts.
    3. Extract CDS features.
    4. Generate CDS coordinates.
    5. Detect overlapping genes.
    6. Classify overlap patterns.

Input:
    GFF3 genome annotations.

Outputs:
    Longest transcript table.
    CDS coordinate table.
    Overlapping gene pairs.
    Overlap classification statistics.
"""

import os
from extract_gff_information import save_to_txt as save_gff_infor
from extract_longest_transcript import process_files as process_transcript
from extract_cds import parse_attributes as extract_cds
from identify_gene_coordinates_via_longest_transcript_cds import parse_attributes as find_longest_cds
from identify_overlapping_genes import process_files as process_overlaps


def main(input_dir, output_dir):
    """
    Run the overlapping-gene identification workflow for all animal genome annotations.
    """

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


            longest_tra_result_file = os.path.join(result_folder, 'protein_gene_info.txt_longest_tra.csv')
            protein_gene_info_file = os.path.join(result_folder, 'protein_gene_info.txt')


            print(f'Processing {gff_file_path}...')

            # Step 1. Extract transcript annotations
            save_gff_infor(gff_file_path , result_folder)
            print(f'Information is saved to {os.path.join(result_folder, "protein_gene_info.txt")}')

            # Step 2. Select representative transcripts
            process_transcript(protein_gene_info_file, result_folder)

            # Step 3. Extract CDS features
            extract_cds(protein_gene_info_file, longest_tra_result_file, '/public2/home/liulm/plant/result1.txt')

            # Step 4. Generate CDS coordinates
            find_longest_cds('/public2/home/liulm/plant/result1.txt', longest_tra_result_file, '/public2/home/liulm/plant/result2.txt')
            process_overlaps('/public2/home/liulm/plant/result2.txt', result_folder)

            print(f'Finished processing {gff_file_path}\n')


if __name__ == "__main__":

    input_dir = '/public2/home/liulm/plant/plant_1'
    output_dir = '/public2/home/liulm/plant/result1_1/'


    main(input_dir, output_dir)
