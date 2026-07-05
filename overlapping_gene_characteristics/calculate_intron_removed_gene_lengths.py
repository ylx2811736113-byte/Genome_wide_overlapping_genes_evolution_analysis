import os
import pandas as pd

# Define the root directory
base_dir = "D:/Work/overlapping_gene/results/animal_new"  # Replace with your actual root directory


def process_file(file_path, output_column):
    """
    Calculate intron-removed gene lengths for a single input file.
    """

    # Read the input table
    df = pd.read_csv(file_path)

    # Check whether the required columns are present
    required_cols = {"GeneID", "gene_length", "intron_length"}

    if not required_cols.issubset(df.columns):
        missing = required_cols - set(df.columns)
        raise ValueError(f"Missing required columns: {missing}")

    # Calculate the total intron length for each gene
    intron_sum = (
        df.groupby("GeneID")["intron_length"]
        .sum()
        .reset_index()
        .rename(columns={"intron_length": "total_intron_length"})
    )

    # Retain one representative record for each gene
    first_records = df.groupby("GeneID", as_index=False).first()

    # Merge gene information with total intron lengths
    merged = pd.merge(
        first_records,
        intron_sum,
        on="GeneID",
        how="inner",
    )

    # Calculate intron-removed gene lengths
    merged[output_column] = (
        merged["gene_length"]
        - merged["total_intron_length"]
    )

    # Remove unnecessary columns
    merged.drop(
        columns=["intron_ID", "intron_length"],
        inplace=True,
        errors="ignore",
    )

    return merged


def main():
    """
    Process all species directories.
    """

    # Traverse all immediate subdirectories
    for folder_name in os.listdir(base_dir):

        folder_path = os.path.join(base_dir, folder_name)

        # Process directories only
        if not os.path.isdir(folder_path):
            continue

        og_file = os.path.join(
            folder_path,
            "og_matched_genes_intron.csv",
        )

        nonog_file = os.path.join(
            folder_path,
            "og_unmatched_genes_intron.csv",
        )

        # Skip if either input file is missing
        if not (
            os.path.exists(og_file)
            and os.path.exists(nonog_file)
        ):
            print(f"Missing input files in {folder_name}")
            continue

        try:
            # Calculate intron-removed gene lengths
            og_result = process_file(
                og_file,
                "OG_Gene_Length_Without_Introns",
            )

            nonog_result = process_file(
                nonog_file,
                "NonOG_Gene_Length_Without_Introns",
            )

            # Combine the results
            combined = pd.concat(
                [og_result, nonog_result],
                ignore_index=True,
            )

            # Save the output table
            output_file = os.path.join(
                folder_path,
                "genes_remove_intron_length.csv",
            )

            combined.to_csv(
                output_file,
                index=False,
                encoding="utf-8",
            )

            print(f"Processed: {folder_name}")

        except Exception as e:
            print(f"Error processing {folder_name}: {e}")


if __name__ == "__main__":
    main()
