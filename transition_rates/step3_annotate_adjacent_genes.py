"""
Step 3. Annotate genes with adjacency information using species-specific intergenic distance thresholds.

Input:
    1. Species-specific ortholog tables.
    2. Gene coordinate tables.
    3. Distance threshold configuration file.

Output:
    Updated ortholog tables containing an additional
    'Position' column indicating whether each gene has
    an adjacent neighboring gene within the specified
    distance threshold.
"""

import os
import pandas as pd


# Root directory containing species-specific results
root_dir = 'D:/Work/overlapping_gene/results/animal_new'

# Distance-threshold configuration table
config_file = (
    'D:/Work/overlapping_gene/results/animal_new/intergenic_distance_thresholds.csv'
)

if not os.path.exists(config_file):
    raise FileNotFoundError(
        f"Configuration file not found:\n{config_file}"
    )

config_df = pd.read_csv(config_file)

required_columns = ['Folder_Name', '5Value']

missing_columns = [
    col for col in required_columns
    if col not in config_df.columns
]

if missing_columns:
    raise ValueError(
        "Missing required columns in configuration file: "
        + ", ".join(missing_columns)
    )

# Build species-to-threshold mapping
threshold_dict = dict(
    zip(
        config_df['Folder_Name'],
        config_df['5Value']
    )
)

processed_species = 0
processed_files = 0
skipped_files = 0

# Process each species
for species in os.listdir(root_dir):

    species_dir = os.path.join(
        root_dir,
        species
    )

    if not os.path.isdir(species_dir):
        continue

    if species not in threshold_dict:

        print(
            f"Skipping {species}: "
            f"distance threshold not found."
        )

        skipped_files += 1
        continue

    distance_threshold = threshold_dict[species]

    print(
        f"\nProcessing {species}"
    )

    print(
        f"Distance threshold = "
        f"{distance_threshold}"
    )

    processed_species += 1

    for subdir, _, files in os.walk(species_dir):

        orthologs_file = next(
            (
                os.path.join(subdir, f)
                for f in files
                if f.endswith(
                    f'orthologues_{species}_final_unique.tsv'
                )
            ),
            None
        )

        genomic_file = next(
            (
                os.path.join(subdir, f)
                for f in files
                if f.endswith(
                    'extracted_genes.csv'
                )
            ),
            None
        )

        if orthologs_file is None or genomic_file is None:
            continue

        try:

            ortholog_df = pd.read_csv(
                orthologs_file,
                sep='\t',
                header=None,
                names=[
                    'c1',
                    'c2',
                    'Protein1',
                    'c3',
                    'Protein2',
                    'c4'
                ]
            )

            genomic_df = pd.read_csv(
                genomic_file
            )

        except Exception as e:

            print(
                f"Skipping {subdir}: "
                f"{e}"
            )

            skipped_files += 1
            continue

        required_gene_columns = [
            'cds_id',
            'Start',
            'End',
            'chr_id'
        ]

        missing = [
            c for c in required_gene_columns
            if c not in genomic_df.columns
        ]

        if missing:

            print(
                f"Skipping {genomic_file}: "
                f"missing columns "
                f"{', '.join(missing)}"
            )

            skipped_files += 1
            continue

        ortholog_df['Position'] = 'not adjacent'

        # Determine whether each gene has neighboring genes within the specified distance threshold
        for index, row in ortholog_df.iterrows():

            protein_id = row['Protein2']

            match = genomic_df[
                genomic_df['cds_id'] == protein_id
            ]

            if match.empty:
                continue

            start = match.iloc[0]['Start']
            end = match.iloc[0]['End']
            chr_id = match.iloc[0]['chr_id']

            previous_genes = genomic_df[
                (genomic_df['chr_id'] == chr_id) &
                (genomic_df['End'] < start) &
                (
                    genomic_df['End']
                    >= start - distance_threshold
                )
            ]

            next_genes = genomic_df[
                (genomic_df['chr_id'] == chr_id) &
                (genomic_df['Start'] > end) &
                (
                    genomic_df['Start']
                    <= end + distance_threshold
                )
            ]

            if (
                not previous_genes.empty
                or
                not next_genes.empty
            ):

                ortholog_df.at[
                    index,
                    'Position'
                ] = 'adjacent'

        output_file = os.path.join(
            subdir,
            f'protein_with_position_update_5Value_{distance_threshold}.csv'
        )

        ortholog_df.to_csv(
            output_file,
            index=False
        )

        processed_files += 1

        print(
            f"Saved: {output_file}"
        )

print("\nFinished.")

print(
    f"Species processed : {processed_species}"
)

print(
    f"Files processed   : {processed_files}"
)

print(
    f"Files skipped     : {skipped_files}"
)
