import os
import pandas as pd

# Define the root directory and output file path
root_dir = "D:/Work/overlapping_gene/results/animal_new"  # Replace with your actual root directory
output_file = os.path.join(
    root_dir,
    "summary_gene_and_intron_lengths.csv"
)


def calculate_mean(file_path, column_name):
    """
    Calculate the mean value of a specified column in a CSV file.
    """
    if not os.path.exists(file_path):
        return float("nan")

    df = pd.read_csv(file_path)

    if column_name not in df.columns:
        print(f"Missing column '{column_name}' in {file_path}")
        return float("nan")

    return df[column_name].mean()


def summarize_species(folder_path):
    """
    Summarize average gene length and intron length for one species.
    """

    return {
        "Species": os.path.basename(folder_path),

        "Average_OG_Intron_Length":
            calculate_mean(
                os.path.join(folder_path, "og_matched_genes_intron.csv"),
                "intron_length"
            ),

        "Average_NonOG_Intron_Length":
            calculate_mean(
                os.path.join(folder_path, "og_unmatched_genes_intron.csv"),
                "intron_length"
            ),

        "Average_OG_Gene_Length":
            calculate_mean(
                os.path.join(folder_path, "og_matched_genes.csv"),
                "gene_length"
            ),

        "Average_NonOG_Gene_Length":
            calculate_mean(
                os.path.join(folder_path, "og_unmatched_genes.csv"),
                "gene_length"
            ),
    }


def main():
    """
    Traverse all species directories and summarize the statistics.
    """

    results = []

    required_files = {
        "og_matched_genes_intron.csv",
        "og_unmatched_genes_intron.csv",
        "og_matched_genes.csv",
        "og_unmatched_genes.csv",
    }

    # Traverse all subdirectories
    for subdir, _, files in os.walk(root_dir):

        # Process directories containing all required files
        if required_files.issubset(set(files)):

            print(f"Processing: {subdir}")

            try:
                results.append(
                    summarize_species(subdir)
                )

            except Exception as e:
                print(f"Error processing {subdir}: {e}")

    # Save the summary table
    result_df = pd.DataFrame(results)

    result_df.to_csv(
        output_file,
        index=False,
        encoding="utf-8"
    )

    print(f"Results have been saved to {output_file}")


if __name__ == "__main__":
    main()
