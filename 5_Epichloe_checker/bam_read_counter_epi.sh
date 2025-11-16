#!/bin/bash

# Folder containing your BAM files
bam_folder="/scratch/drt83172/Wallace_lab/RNA_SEQ/filtered_bams_epi/sep"

# Go to the folder
cd "$bam_folder" || exit

# Loop over all BAM files and print filename + read count
for bam in *.bam; do
    count=$(samtools view -c "$bam")
    echo -e "$bam\t$count"
done > /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/Lists/bam_read_counts.txt