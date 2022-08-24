#!/bin/bash
# Start by activating decona

# Usage: bash decona.sh <path.to.out.folder>

OUTPUT="${1}"

# Grab the parameters file so you have the clustering threshold

source "${OUTPUT}"/banzai_params.sh

# All files of interest should be in noprimers/file/plate/*.fastq
cd "${OUTPUT}"

for folder in noprimers/*; do
  cd ${folder}
  decona -c "${CLUSTER_SIM}" -n 5 -M
  cd "${OUTPUT}"
done
