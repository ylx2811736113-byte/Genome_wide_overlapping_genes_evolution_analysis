"""
Step 12. Summarize transition-rate estimates from BayesTraits output files.

This script extracts transition-rate parameters and
root-state probabilities from the final two lines of
each BayesTraits result file and combines them into
a single summary table.

Input:
    *.txt (BayesTraits output files)

Output:
    transition_rates_animal.csv
"""

import os
import pandas as pd


# Directory containing BayesTraits result files
input_dir = (
    "D:/Work/overlapping_gene/revision/results/"
    "traits_result_animal_013026"
)

# Output summary table
output_file = (
    "D:/Work/overlapping_gene/revision/results/"
    "transition_rates_animal.csv"
)

# Parameters to extract
target_columns = [
    "q01", "q02", "q03",
    "q10", "q12", "q13",
    "q20", "q21", "q23",
    "q30", "q31", "q32",
    "Root P(0)",
    "Root P(1)",
    "Root P(2)",
    "Root P(3)"
]

results = []

processed_files = 0
skipped_files = 0


if not os.path.isdir(input_dir):
    raise FileNotFoundError(
        f"Input directory not found: {input_dir}"
    )


for file_name in sorted(os.listdir(input_dir)):

    if not file_name.endswith(".txt"):
        continue

    file_path = os.path.join(
        input_dir,
        file_name
    )

    try:

        with open(file_path, "r") as f:
            lines = f.readlines()

    except Exception as e:

        print(
            f"Skipping {file_name}: "
            f"{e}"
        )

        skipped_files += 1
        continue

    if len(lines) < 2:

        print(
            f"Skipping {file_name}: "
            f"file contains fewer than two lines."
        )

        skipped_files += 1
        continue

    header = lines[-2].rstrip().split("\t")
    values = lines[-1].rstrip().split("\t")

    if len(header) != len(values):

        print(
            f"Skipping {file_name}: "
            f"header/value length mismatch."
        )

        skipped_files += 1
        continue

    parameter_dict = dict(
        zip(header, values)
    )

    record = {
        "File_Name": file_name
    }

    for parameter in target_columns:
        record[parameter] = parameter_dict.get(
            parameter,
            None
        )

    results.append(record)

    processed_files += 1

    print(
        f"Processed: {file_name}"
    )


summary_df = pd.DataFrame(results)

summary_df.to_csv(
    output_file,
    index=False
)

print("\nFinished.")

print(
    f"Files processed : {processed_files}"
)

print(
    f"Files skipped   : {skipped_files}"
)

print(
    f"Summary saved to: {output_file}"
)
