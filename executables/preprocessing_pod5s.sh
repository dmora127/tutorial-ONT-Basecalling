#!/bin/bash
# Preprocess POD5 files by generating summaries and splitting into subsets based on channel
# Usage: ./preprocessing_pod5s.sh <path_to_pod5_dir>

set -euo pipefail
POD5_DIR="$1"

# Generate summary of POD5 files
pod5 view $POD5_DIR --include "read_id, channel" --output summary.tsv
echo "Generated POD5 summary."

# Split POD5 files into subsets based on channel
pod5 subset $POD5_DIR --summary summary.tsv --columns channel --output split_pod5_subsets
echo "Completed splitting POD5 files into subsets based on channel."
echo "Preprocessing of POD5 files completed."

ls split_pod5_subsets > list_of_pod5_files.txt
echo "List of split POD5 files saved to list_of_pod5_files.txt."