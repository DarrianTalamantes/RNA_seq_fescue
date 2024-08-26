#!/bin/bash
# Purpose: This script is here because I want to run FASTQC on things that are not part of the snakemake pipline
# Author: Darrian Talamantes (drt83172@uga.edu)


# Check if a directory is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <directory>"
  exit 1
fi

# Directory containing the files
DIR=$1

# Output directory for FastQC results
OUTPUT_DIR="${DIR}/fastqc_results"
mkdir -p $OUTPUT_DIR

# Run FastQC on all fastq files in the specified directory
>commands.txt
for file in "$DIR"/*.fq.gz; do
  if [ -f "$file" ]; then
    echo "Running FastQC on $file"
    echo "fastqc -o $OUTPUT_DIR $file" >> commands.txt 
  else
    echo "No FASTQ files found in $DIR"
    exit 1
  fi
done

parallel --jobs 12 < commands.txt
echo "FastQC analysis complete. Results are in $OUTPUT_DIR"