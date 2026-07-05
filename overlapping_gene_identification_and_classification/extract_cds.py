"""
Extract CDS annotations for the longest transcript of each protein-coding gene.

This script parses a genome annotation (GFF3) file and retrieves all CDS
features associated with the longest transcript of each protein-coding gene.
The resulting CDS coordinates and corresponding protein identifiers are
exported as an intermediate table for downstream overlap analyses.

Input:
    1. Genome annotation file (GFF3)
    2. Longest transcript table

Output:
    CDS annotation table containing chromosome, genomic coordinates,
    strand, gene ID, transcript ID, CDS ID, and protein ID.
"""

import pandas as pd
import re


def parse_attributes(gff_file, trafile_path, resultfile_path):
    """
    Extract CDS annotations associated with the longest transcript
    of each protein-coding gene and export the results.
    """

    # Load transcript information
    tra_df = pd.read_csv(trafile_path, sep='\t', header=0)

    # Read animal GFF3 file
    gff_df = pd.read_csv(
        gff_file,
        sep='\t',
        comment='#',
        header=None,
        names=[
            'chr_id', 'seqid', 'source', 'type',
            'start', 'end', 'score', 'strand',
            'phase', 'attributes'
        ]
    )

    # Alternative parser for plant GFF3 files
    # gff_df = pd.read_csv(
    #     gff_file,
    #     sep='\t',
    #     comment='#',
    #     header=None,
    #     names=[
    #         'chr_id', 'source', 'type',
    #         'start', 'end', 'score',
    #         'strand', 'phase', 'attributes'
    #     ]
    # )

    # Initialize result DataFrame
    result_df = pd.DataFrame(
        columns=[
            'chr_id',
            'type',
            'start',
            'end',
            'strand',
            'gene_id',
            'transcript_id',
            'cds_id',
            'protein_id'
        ]
    )

    # Iterate over transcripts
    for i in range(len(tra_df)):

        tra_id = tra_df.loc[i, 'transcript_ID']
        gene_id = tra_df.loc[i, 'GeneID']

        # Retrieve CDS features associated with the transcript
        matching_rows = gff_df[
            (gff_df['type'] == 'CDS') &
            (gff_df['attributes'].str.contains(
                f'Parent=rna-{tra_id}',
                na=False
            ))
        ]

        # Extract CDS and protein identifiers
        for _, row in matching_rows.iterrows():

            cds_id_match = re.search(
                r'ID=([^;]+);?',
                row['attributes']
            )

            # Animal GFF3 format
            cds_id = (
                cds_id_match.group(1).lstrip('cds-')
                if cds_id_match else None
            )

            # Plant GFF3 format
            # cds_id = (
            #     cds_id_match.group(1)
            #     if cds_id_match else None
            # )

            protein_match = re.search(
                r'protein_id=([^;]+);?',
                row['attributes']
            )

            protein_id = (
                protein_match.group(1)
                if protein_match else None
            )

            result_row = {
                'chr_id': row['chr_id'],
                'type': row['type'],
                'start': row['start'],
                'end': row['end'],
                'strand': row['strand'],
                'gene_id': gene_id,
                'transcript_id': tra_id,
                'cds_id': cds_id,
                'protein_id': protein_id
            }

            result_df = pd.concat(
                [result_df, pd.DataFrame([result_row])],
                ignore_index=True
            )

            # Print extracted record
            print(result_row)

    # Save results
    result_df.to_csv(
        resultfile_path,
        sep=' ',
        index=False
    )


# Clear previous output file
with open('/public2/home/liulm/plant/result1.txt', 'a') as f1:
    f1.truncate(0)
