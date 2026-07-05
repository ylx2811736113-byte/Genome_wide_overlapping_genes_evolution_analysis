#!/bin/bash

# Define the tree file and the directory containing the trait files
TREE_FILE="/public2/home/liulm/overlapping_gene/plant_020426.trees"
TRAITS_DIR="/public2/home/liulm/overlapping_gene/traits_files_plant_020426"
OUTPUT_DIR="/public2/home/liulm/overlapping_gene/traits_result_plant_020426"

# Create the output directory
mkdir -p $OUTPUT_DIR

# Iterate through all trait files
for TRAIT_FILE in $TRAITS_DIR/*_traits.txt
do
    # Extract the trait file name without the directory path
    BASENAME=$(basename $TRAIT_FILE)

    # Set the log file name with the .Log.txt extension
    LOG_FILE="${TRAIT_FILE}.Log.txt"

    # Run BayesTraits analysis by piping the interactive commands into the program
    echo -e "1\n1\nrun\n" | /public2/home/liulm/BayesTraitsV4.1.2-Linux/BayesTraitsV4 $TREE_FILE $TRAIT_FILE

    # Check whether the log file has been generated
    if [ -f "$LOG_FILE" ]; then
        # Copy the log file to the output directory
        cp "$LOG_FILE" "$OUTPUT_DIR/"
        echo "Processed $TRAIT_FILE, log saved as $LOG_FILE and copied to $OUTPUT_DIR"
    else
        echo "Log file for $TRAIT_FILE was not found."
    fi
done
