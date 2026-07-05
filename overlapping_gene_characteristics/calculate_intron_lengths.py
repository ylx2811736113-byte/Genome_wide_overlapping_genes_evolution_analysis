import os
import pandas as pd

# Define the root directory
root_dir = "D:/Work/overlapping_gene/results/animal_new"  # Replace with your actual root directory


def process_file(file_path):
    """
    Calculate intron lengths from exon coordinates for each transcript.
    """

    # Read the exon annotation table
    df = pd.read_csv(file_path, sep="\t")

    # Check whether the required columns are present
    required_cols = {"transcript_ID", "start", "end"}

    if not required_cols.issubset(df.columns):
        missing = required_cols - set(df.columns)
        raise ValueError(f"Missing required columns: {missing}")

    intron_records = []

    # Process each transcript independently
    for transcript_id, group in df.groupby("transcript_ID"):

        # Sort exons by genomic coordinates
        group = (
            group.sort_values(["start", "end"])
            .reset_index(drop=True)
        )

        # Calculate intron lengths between adjacent exons
        for i in range(1, len(group)):

            intron_length = (
                group.loc[i, "start"]
                - group.loc[i - 1, "end"]
                - 1
            )

            intron_records.append(
                {
                    "transcript_ID": transcript_id,
                    "intron_ID": f"intron{i}",
                    "intron_length": intron_length,
                }
            )

    return pd.DataFrame(intron_records)


def main():
    """
    Traverse all species directories and calculate intron lengths.
    """

    # Traverse all subdirectories
    for subdir, _, _ in os.walk(root_dir):

        input_file = os.path.join(subdir, "exon.txt")

        if not os.path.exists(input_file):
            continue

        print(f"Processing: {input_file}")

        try:

            result = process_file(input_file)

            output_file = os.path.join(
                subdir,
                "intron.txt",
            )

            result.to_csv(
                output_file,
                sep="\t",
                index=False,
                encoding="utf-8",
            )

            print(f"Saved: {output_file}")

        except Exception as e:

            print(f"Error processing {input_file}: {e}")


if __name__ == "__main__":
    main()
