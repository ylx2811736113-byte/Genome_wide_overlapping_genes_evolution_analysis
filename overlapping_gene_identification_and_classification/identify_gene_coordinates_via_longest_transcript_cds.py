"""
Generate CDS-level coordinates for representative transcripts.

Input:
    1. CDS annotation table
    2. Representative transcript table

Output:
    CDS coordinate table for downstream overlapping-gene detection.
"""

import pandas as pd


def parse_attributes(cds_file_path, trafile_path, resultfile_path):
    """
    Generate CDS coordinates for the representative transcript of each gene.
    """

    # Load CDS and transcript annotations
    cds_df = pd.read_csv(cds_file_path, sep=' ', header=0)
    tra_df = pd.read_csv(trafile_path, sep='\t', header=0)

    # Output table
    result_df = pd.DataFrame(columns=['chr_id', 'Type', 'Start', 'End', 'Strand', 'GeneID', 'transcript_id', 'cds_id', 'Protein_id'])
    #result_df = pd.DataFrame(
        #columns=['chr_id', 'Type', 'Start', 'End', 'Strand', 'GeneID', 'transcript_id', 'cds_id'])

    # Process each transcript
    for i in range(len(tra_df)):
        tra_id = tra_df.loc[i, 'transcript_ID']
        gene_id = tra_df.loc[i, 'GeneID']

        # Retrieve CDS fragments belonging to the transcript
        matching_rows = cds_df[cds_df['transcript_id'] == tra_id]

        if not matching_rows.empty:
            # Define CDS boundaries
            start = matching_rows['start'].min()
            end = matching_rows['end'].max()

            example_row = matching_rows.iloc[0]

            result_row = {
                'chr_id': example_row['chr_id'],
                'Type': 'CDS',
                'Start': start,
                'End': end,
                'Strand': example_row['strand'],
                'GeneID': gene_id,
                'transcript_id': tra_id,
                'cds_id': example_row['cds_id'],
                'Protein_id': example_row['protein_id']
            }
            result_df = pd.concat([result_df, pd.DataFrame([result_row])], ignore_index=True)

    # Export results
    result_df.to_csv(resultfile_path, sep=' ', index=False)


with open('/public2/home/liulm/plant/result2.txt', 'a') as f1:
    f1.truncate(0)



