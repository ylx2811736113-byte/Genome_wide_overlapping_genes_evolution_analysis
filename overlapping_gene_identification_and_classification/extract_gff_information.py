"""
Extract protein-coding gene annotations from an animal GFF/GFF3 file.

This script parses a genome annotation (GFF/GFF3) file and extracts
protein-coding gene features together with their corresponding
mRNA and CDS records. Chromosome identifiers are retrieved from
the region annotations and appended to each extracted record.
The resulting annotation table is exported for downstream
identification of overlapping genes.

Input:
    Genome annotation file (GFF/GFF3)

Output:
    protein_gene_info.txt containing chromosome ID, genomic
    coordinates, strand, feature type, and annotation attributes
    for protein-coding genes, mRNAs, and CDSs.
"""

import os


def extract_gene_mrna_cds_info_with_chr(gff_file):
    """
    Extract protein-coding gene, mRNA, and CDS records from a GFF/GFF3 file and associate them with chromosome identifiers.
    """

    gene_info = []
    current_chr_id = None
    capture = False

    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')
            if len(columns) < 9:
                continue

            seq_id, source, feature_type, start, end, score, strand, phase, attributes = columns

            if feature_type == 'region':
                capture = False  # Reset capture flag for each new region
                # Extract chromosome id from attributes
                for attribute in attributes.split(';'):
                    if attribute.startswith('Name='):
                        current_chr_id = attribute.split('=')[1]
                        break

            elif feature_type == 'gene' and 'protein_coding' in attributes:
                if current_chr_id:
                    gene_info.append([current_chr_id] + columns)
                capture = True  # Start capturing mRNA and CDS

            elif capture and feature_type in ['mRNA', 'CDS']:
                if current_chr_id:
                    gene_info.append([current_chr_id] + columns)

            # If we encounter a new region or gene, stop capturing unless it meets criteria
            if feature_type == 'region' or feature_type =='pseudogene' or (feature_type == 'gene' and 'protein_coding' not in attributes):
                capture = False

    return gene_info


def save_to_txt(gff_f, output_dir):
    """
    Save the extracted protein-coding annotation records to a tab-delimited text file.
    """

    gene_mrna_cds_info = extract_gene_mrna_cds_info_with_chr(gff_f)
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'protein_gene_info.txt')

    with open(output_file, 'w') as file:
        # Write header
        file.write('chr_id\tseq_id\tsource\tfeature_type\tstart\tend\tscore\tstrand\tphase\tattributes\n')
        for entry in gene_mrna_cds_info:
            file.write('\t'.join(entry) + '\n')
