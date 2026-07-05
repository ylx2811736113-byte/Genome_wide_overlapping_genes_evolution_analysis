"""
Extract protein-coding CDS annotations and generate representative CDS coordinates.

Input:
    Genome annotation file in GFF/GFF3 format.

Output:
    1. All protein-coding CDS records.
    2. Representative CDS table containing the longest CDS for each gene.
"""

import csv
import re
import pandas as pd
import os


resultfile_path = '/public2/home/liulm/plant/result.txt'


def parse_attributes(gff_file):
    """
    Extract protein-coding CDS annotations from a GFF3 file.
    """

    gene_biotype_dict = {}

    with open(gff_file, 'r') as file:

        reader = csv.reader(file, delimiter='\t')

        for row in reader:

            if re.match('##', row[0]):
                continue
            elif re.match('#!' or '###', row[0]):
                continue
            elif re.match('#', row[0]):
                continue

            # Store gene biotypes
            if row[2] == 'gene':

                gene_id = re.search(
                    r'ID=([^;]+);?',
                    row[8]
                ).group(1)

                gene_biotype_match = re.search(
                    r'gene_biotype=([^;]+);?',
                    row[8]
                )

                if gene_biotype_match:
                    gene_biotype_dict[gene_id] = (
                        gene_biotype_match.group(1)
                    )

            # Extract protein-coding CDS features
            if row[2] == 'CDS':

                parent_id = re.search(
                    r'Parent=([^;]+);?',
                    row[8]
                ).group(1)

                gene_biotype = gene_biotype_dict.get(
                    parent_id,
                    None
                )

                if gene_biotype == 'protein_coding':

                    gene_id = parent_id.lstrip('gene-')

                    cds_id = re.search(
                        r'ID=([^;]+);?',
                        row[8]
                    ).group(1)

                    cds_id = cds_id.lstrip('cds-')

                    cds_length = (
                        int(row[4]) - int(row[3])
                    )

                    protein_id = re.search(
                        r'protein_id=([^;]+);?',
                        row[8]
                    ).group(1)

                    with open(
                        resultfile_path,
                        'a'
                    ) as resultfile:

                        print(
                            row[0],
                            row[2],
                            row[3],
                            row[4],
                            row[6],
                            gene_id,
                            cds_id,
                            cds_length,
                            protein_id,
                            file=resultfile
                        )

    gene_info_df = pd.read_csv(
        resultfile_path,
        sep=' ',
        header=None,
        names=[
            "Seq_id",
            "Type",
            "Start",
            "End",
            "Strand",
            "GeneID",
            "CDS_ID",
            "CDS_length",
            "Protein_id"
        ],
        usecols=range(9)
    )

    return gene_info_df


def sort_cds_records(genes_info_df):
    """
    Sort CDS records according to genomic coordinates.
    """

    genes_info_df.sort_values(
        by=["Seq_id", "Start"],
        inplace=True
    )

    genes_info_df = (
        genes_info_df
        .reset_index(drop=True)
    )

    genes_info_df['Index'] = (
        genes_info_df.index + 1
    )

    return genes_info_df


def get_longest_cds(group):
    """
    Select the longest CDS for each gene.
    """

    longest = group.loc[
        group['CDS_length'].idxmax()
    ]

    return pd.Series({
        'cds_no': len(group),
        'chr_id': longest['Seq_id'],
        'GeneID': longest['GeneID'],
        'CDS_ID': longest['CDS_ID'],
        'Start': longest['Start'],
        'End': longest['End'],
        'Strand': longest['Strand'],
        'Protein_id': longest['Protein_id']
    })


def process_files(file_paths, output_dir):
    """
    Generate CDS tables for downstream overlapping-gene identification.
    """

    with open(
        '/public2/home/liulm/plant/result.txt',
        'a'
    ) as f1:
        f1.truncate(0)

    cds_df = parse_attributes(file_paths)

    sorted_df = sort_cds_records(cds_df)

    representative_cds = (
        sorted_df
        .groupby('GeneID')
        .apply(get_longest_cds)
        .reset_index(drop=True)
    )

    representative_cds.sort_values(
        by=['chr_id', 'Start'],
        inplace=True
    )

    extracted_cds_file = os.path.join(
        output_dir,
        os.path.basename(file_paths)
        + '_extracted_cds.csv'
    )

    longest_cds_file = os.path.join(
        output_dir,
        os.path.basename(file_paths)
        + '_longest_cds.csv'
    )

    sorted_df.to_csv(
        extracted_cds_file,
        sep='\t',
        index=False
    )

    representative_cds.to_csv(
        longest_cds_file,
        sep='\t',
        index=False
    )

    print(
        f"Extracted CDS records saved to "
        f"{extracted_cds_file}"
    )

    print(
        f"Representative CDS table saved to "
        f"{longest_cds_file}"
    )


# Reset temporary file
with open(
    '/public2/home/liulm/plant/result.txt',
    'a'
) as f1:
    f1.truncate(0)