#!/bin/bash
# Run Dorado on the EP for each POD5 file (non-resumeable)
# Usage: ./run_dorado.sh "<dorado_args>" <dorado_model_name> <input_pod5_file>

set -euo pipefail

DORADO_ARGS="$1"
DORADO_MODEL_NAME="$2"
INPUT_POD5="$3"
BAM_FILE="${INPUT_POD5}.bam"
FASTQ_FILE="${BAM_FILE}.fastq"

echo "Running Dorado with args: ${DORADO_ARGS}"
echo "Model Name: ${DORADO_MODEL_NAME}"
echo "Input POD5: ${INPUT_POD5}"
echo "Output BAM: ${BAM_FILE}"

# untar your Dorado basecalling models
tar -xvzf ${DORADO_MODEL_NAME}.tar.gz
rm ${DORADO_MODEL_NAME}.tar.gz
echo "Completed model extraction."

# Run Dorado
dorado ${DORADO_ARGS} "${DORADO_MODEL_NAME}" "${INPUT_POD5}" > "${BAM_FILE}"
echo "Completed Dorado basecalling."

# Convert BAM to FASTQ
samtools fastq "${BAM_FILE}" > "${FASTQ_FILE}"
echo "Completed BAM → FASTQ conversion."
