"""
Identify and classify overlapping genes.

Input:
    CDS coordinate table of representative transcripts.

Output:
    1. Sorted gene coordinate table.
    2. Overlapping gene pairs.
    3. Summary statistics of overlap categories.

Overlap categories:
    General overlap types:
        All++, All--, All+-, All-+
        Part++, Part--, Part+-, Part-+

    Narrow overlap types:
        Nested
        Co-oriented
        Convergent
        Divergent
"""

import pandas as pd
import os



def parse_attributes(cds_info_path):
    """
    Load CDS coordinate information.
    """
    gene_info_df = pd.read_csv(cds_info_path, sep=' ', header=0)
    return gene_info_df


def extract_gene_info(genes_info_df):
    """
    Sort genes according to genomic coordinates and assign indices.
    """
    genes_info_df.sort_values(by=["chr_id", "Start"], inplace=True)
    # genes_info_df['Attributes_Dict'] = genes_info_df['Attributes'].apply(parse_attributes)
    # genes_info_df['Gene_biotype'] = genes_info_df.apply(lambda x: x['Attributes_Dict'].get('locus_type', 'protein_coding'), axis=1)

    genes_info_df = genes_info_df.reset_index(drop=True)
    genes_info_df['Index'] = genes_info_df.index + 1
    return genes_info_df


def identify_and_describe_overlaps(df):
    """
    Detect overlapping gene pairs and classify overlap patterns.
    """
    df = df.sort_values(by=["chr_id", "Start"]).reset_index(drop=True)
    df['Index'] = df.index + 1
    overlaps_info = []
    overlap_count = 0
    goverlap_type = 'Not determined'
    noverlap_type = 'Not determined'
    for i in range(len(df) - 1):
        for j in range(i + 1, len(df)):
            # Candidate overlap on the same chromosome
            if df.loc[i, 'chr_id'] == df.loc[j, 'chr_id'] and df.loc[j, 'Start'] <= df.loc[i, 'End']:
                overlap_count += 1
                # Same-strand overlaps (+/+)
                if df.loc[i, 'Strand'] == '+' and df.loc[j, 'Strand'] == '+':
                    if df.loc[i, 'End'] >= df.loc[j, 'End']:
                        goverlap_type = 'All++'
                        noverlap_type = 'Nested'
                    elif df.loc[i, 'Start'] < df.loc[j, 'Start'] and df.loc[i, 'End'] < df.loc[j, 'End']:
                        goverlap_type = 'Part++'
                        noverlap_type = 'Co-oriented'
                    else:
                        goverlap_type = 'All++'
                        noverlap_type = 'Nested'
                # Same-strand overlaps (-/-)
                if df.loc[i, 'Strand'] == '-' and df.loc[j, 'Strand'] == '-':
                    if df.loc[i, 'End'] >= df.loc[j, 'End']:
                        goverlap_type = 'All--'
                        noverlap_type = 'Nested'
                    elif df.loc[i, 'Start'] < df.loc[j, 'Start'] and df.loc[i, 'End'] < df.loc[j, 'End']:
                        goverlap_type = 'Part--'
                        noverlap_type = 'Co-oriented'
                    else:
                        goverlap_type = 'All--'
                        noverlap_type = 'Nested'
                # Opposite-strand overlaps (+/-)
                if df.loc[i, 'Strand'] == '+' and df.loc[j, 'Strand'] == '-':
                    if df.loc[i, 'End'] >= df.loc[j, 'End']:
                        goverlap_type = 'All+-'
                        noverlap_type = 'Nested'
                    elif df.loc[i, 'Start'] < df.loc[j, 'Start'] and df.loc[i, 'End'] < df.loc[j, 'End']:
                        goverlap_type = 'Part+-'
                        noverlap_type = 'Convergent'
                    else:
                        goverlap_type = 'All+-'
                        noverlap_type = 'Nested'
                # Opposite-strand overlaps (-/+)
                if df.loc[i, 'Strand'] == '-' and df.loc[j, 'Strand'] == '+':
                    if df.loc[i, 'End'] >= df.loc[j, 'End']:
                        goverlap_type = 'All-+'
                        noverlap_type = 'Nested'
                    elif df.loc[i, 'Start'] < df.loc[j, 'Start'] and df.loc[i, 'End'] < df.loc[j, 'End']:
                        goverlap_type = 'Part-+'
                        noverlap_type = 'Divergent'
                    else:
                        goverlap_type = 'All-+'
                        noverlap_type = 'Nested'
                overlaps_info.append({
                    'Overlap_No': overlap_count,
                    'Gene1_Index': df.loc[i, 'Index'],
                    'Gene1_chr_id': df.loc[i, 'chr_id'],
                    'Gene1_ID': df.loc[i, 'GeneID'],
                    'Gene1_Start': df.loc[i, 'Start'],
                    'Gene1_End': df.loc[i, 'End'],
                    'Gene1_Strand': df.loc[i, 'Strand'],
                    'Gene1_protein_id': df.loc[i, 'Protein_id'],
                    'Gene2_Index': df.loc[j, 'Index'],
                    'Gene2_chr_id': df.loc[j, 'chr_id'],
                    'Gene2_ID': df.loc[j, 'GeneID'],
                    'Gene2_Start': df.loc[j, 'Start'],
                    'Gene2_End': df.loc[j, 'End'],
                    'Gene2_Strand': df.loc[j, 'Strand'],
                    'Gene2_protein_id': df.loc[j, 'Protein_id'],
                    'General_overlap_Type': goverlap_type,
                    'Narrow_overlap_Type': noverlap_type
                })
    return pd.DataFrame(overlaps_info)


def generate_overlaps_statistics(overlap_df, output_file):
    """
    Summarize overlap categories.
    """
    goverlap_types_count = overlap_df['General_overlap_Type'].value_counts()
    noverlap_types_count = overlap_df['Narrow_overlap_Type'].value_counts()

    with open(output_file, 'w') as f:
        f.write("Overlap Categories:\n")
        f.write(goverlap_types_count.to_string())
        f.write(noverlap_types_count.to_string())


def process_files(file_paths, output_dir):
    """
    Run overlap detection and export results.
    """
    geneinfo_df = parse_attributes(file_paths)
    extracted_df = extract_gene_info(geneinfo_df)
    overlap_df = identify_and_describe_overlaps(extracted_df)
    extracted_genes_file = os.path.join(output_dir, os.path.basename(file_paths) + '_extracted_genes.csv')
    overlap_genes_file = os.path.join(output_dir, os.path.basename(file_paths) + '_overlap_genes.csv')
    extracted_df.to_csv(extracted_genes_file, sep='\t', index=False)
    overlap_df.to_csv(overlap_genes_file, sep='\t', index=False)

    statistics_file = os.path.join(output_dir, os.path.basename(file_paths) + '_statistics.txt')
    generate_overlaps_statistics(overlap_df, statistics_file)
    print(f"  Extracted genes saved to {extracted_genes_file}")
    print(f"  Overlap genes saved to {overlap_genes_file}")
    print(f"  Statistics saved to {statistics_file}")












