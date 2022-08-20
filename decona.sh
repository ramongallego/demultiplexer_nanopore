#!/bin/bash
# Start by activating decona

# Usage: bash decona.sh <path.to.out.folder>

OUTPUT_DIRECTORY="${1}"

# All files of interest should be in noprimers/file/plate/*.fastq
cd "${OUTPUT_DIRECTORY}"

for folder in noprimers/*/*; do
  cd ${folder}
  decona -c 0.85 -n 10 -M
  cd "${OUTPUT_DIRECTORY}"
done
