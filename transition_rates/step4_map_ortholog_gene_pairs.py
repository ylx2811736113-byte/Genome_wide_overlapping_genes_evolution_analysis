"""
Step 4. Map reference gene pairs to orthologous genes in target species.

Input:
    1. Reference gene-pair table.
    2. Species-specific ortholog mapping tables.

Output:
    Ortholog gene-pair tables for each species.
"""

import os
import csv
from collections import defaultdict


# Reference gene-pair table
reference_pairs = (
    'D:/Work/overlapping_gene/results/model_species/Homo_sapiens/gene_pairs_with_match.csv'
)

# Root directory containing species-specific results
root_folder = (
    'D:/Work/overlapping_gene/results/animal_new'
)

if not os.path.exists(reference_pairs):
    raise FileNotFoundError(
        f"Reference file not found:\n{reference_pairs}"
    )

processed_files = 0
skipped_files = 0

# Process each ortholog mapping table
for folder_name, _, filenames in os.walk(root_folder):

    for filename in filenames:

        if not filename.endswith('final_unique.tsv'):
            continue

        ortholog_file = os.path.join(
            folder_name,
            filename
        )

        output_file = os.path.join(
            folder_name,
            'ortholog_pairs.csv'
        )

        print(
            f"\nProcessing: {ortholog_file}"
        )

        ortholog_dict = defaultdict(list)

        # Read ortholog mapping table
        try:

            with open(
                ortholog_file,
                'r',
                newline=''
            ) as f:

                reader = csv.reader(
                    f,
                    delimiter='\t'
                )

                for row in reader:

                    if len(row) < 6:
                        continue

                    reference_gene = row[1]
                    target_gene = row[4]

                    if (
                        reference_gene
                        and
                        target_gene
                    ):

                        ortholog_dict[
                            reference_gene
                        ].append(
                            target_gene
                        )

        except Exception as e:

            print(
                f"Skipping {ortholog_file}: "
                f"{e}"
            )

            skipped_files += 1
            continue

        # Generate mapped gene-pair table
        try:

            with open(
                reference_pairs,
                'r',
                newline=''
            ) as infile, open(
                output_file,
                'w',
                newline=''
            ) as outfile:

                reader = csv.reader(infile)
                writer = csv.writer(outfile)

                for row in reader:

                    if len(row) < 2:
                        continue

                    gene1, gene2 = row[:2]

                    ortholog1 = ";".join(
                        ortholog_dict.get(
                            gene1,
                            ['no']
                        )
                    )

                    ortholog2 = ";".join(
                        ortholog_dict.get(
                            gene2,
                            ['no']
                        )
                    )

                    writer.writerow([

                        gene1,
                        ortholog1,

                        gene2,
                        ortholog2

                    ])

        except Exception as e:

            print(
                f"Failed to generate "
                f"{output_file}: {e}"
            )

            skipped_files += 1
            continue

        processed_files += 1

        print(
            f"Saved: {output_file}"
        )

print("\nFinished.")

print(
    f"Files processed : {processed_files}"
)

print(
    f"Files skipped   : {skipped_files}"
)
