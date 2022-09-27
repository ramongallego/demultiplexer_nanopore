#!/bin/bash
# Start by activating decona

# Usage: bash decona.sh <path.to.params> <path.to.folder>
# Folder has all the subfolders with the demultiplexed data

OUTPUT="${2}"

# Grab the parameters file so you have the clustering threshold

source "${1}"

# All files of interest should be in noprimers/file/plate/*.fastq
cd "${OUTPUT}"

for folder in */; do
  cd ${folder}  
  decona -l "${MIN_LENGTH}" -m "${MAX_LENGTH}" -q 10 -c "${CLUSTER_SIM}" -n 10 -M
  cd ..
done
