"""
Identify the longest transcript for each protein-coding gene.

This script parses a genome annotation (GFF3) file, extracts all
mRNA annotations, and calculates the genomic length of each
transcript. For genes with multiple transcript isoforms, the
longest transcript is selected as the representative transcript.
The extracted transcript annotations and the longest-transcript
summary are exported for downstream CDS extraction and overlap
analysis.

Input:
    Genome annotation file (GFF3)

Output:
    1. Extracted mRNA annotation table
    2. Longest-transcript summary table
"""

import csv
import re
import pandas as pd
import os


def parse_attributes(gff_file):
    """Extract mRNA annotations and transcript information from a GFF/GFF3 file."""
    with (open(gff_file, 'r') as file):
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if re.match('##', row[0]):
                continue
            elif re.match('#!' or '###', row[0]):
                continue
            elif re.match('#', row[0]):
                continue
            # plant
            if row[2] == 'mRNA':
            # animal
            #if row[3] == 'mRNA':
                # print(row)
                gene_id = re.search(r'Parent=([^;]+);?', row[9]).group(1)
                gene_id = gene_id.lstrip('gene-')
                transcript_id = re.search(r'ID=([^;]+);?', row[9]).group(1)
                transcript_id = transcript_id.lstrip('rna-')
                # plant
                transcript_length = int(row[4]) - int(row[3]) + 1
                # animal
                # transcript_length = int(row[5]) - int(row[4]) + 1
                with open(resultfile_path, 'a') as resultfile:
                    print(row[0], row[3], row[4], row[5], row[7], gene_id, transcript_id, transcript_length,
                          file=resultfile)
    gene_info_df = pd.read_csv(resultfile_path, sep=' ', header=None,
                               names=["Chr_id", "Type", "Start", "End", "Strand", "GeneID", "TranscriptID",
                                      "Transcript_length"],
                               usecols=range(8))
    return gene_info_df


def sort_tra_info(genes_info_df):
    """Sort transcript annotations by chromosome and genomic coordinate."""
    genes_info_df.sort_values(by=["Chr_id", "Start"], inplace=True)
    # Reset index and add an Index column
    genes_info_df = genes_info_df.reset_index(drop=True)
    genes_info_df['Index'] = genes_info_df.index + 1
    return genes_info_df


def get_longest_transcript(group):
    """
    Identify the longest transcript for each gene.
    """
    longest = group.loc[group['Transcript_length'].idxmax()]
    return pd.Series({
        'transcript_no': len(group),  # Number of transcripts associated with this gene
        'transcript_chr_id': longest['Chr_id'],
        'GeneID': longest['GeneID'],
        'transcript_ID': longest['TranscriptID'],
        'Gene_Start': longest['Start'],
        'Gene_End': longest['End'],
        'Gene_Strand': longest['Strand']
    })


def process_files(file_paths, output_dir):
    """Process a GFF3 file, extract transcript annotations, identify the longest transcript for each gene, and export the results."""
    mrnainfo_df = parse_attributes(file_paths)
    sorted_df = sort_tra_info(mrnainfo_df)
    transcript_info = sorted_df.groupby('GeneID').apply(get_longest_transcript).reset_index(drop=True)

    # Sort transcript_info
    transcript_info.sort_values(by=['transcript_chr_id', 'Gene_Start'], inplace=True)

    extracted_mrna_file = os.path.join(output_dir, os.path.basename(file_paths) + '_extracted_mrna.csv')
    longest_transcript_file = os.path.join(output_dir, os.path.basename(file_paths) + '_longest_tra.csv')

    sorted_df.to_csv(extracted_mrna_file, sep='\t', index=False)
    transcript_info.to_csv(longest_transcript_file, sep='\t', index=False)

    print(f"  Extracted mRNA saved to {extracted_mrna_file}")
    print(f"  Longest transcript saved to {longest_transcript_file}")


with open('/public2/home/liulm/plant/result.txt', 'a') as f1:
    f1.truncate(0)
resultfile_path = '/public2/home/liulm/plant/result.txt'

