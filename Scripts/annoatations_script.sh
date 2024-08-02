#!/bin/bash
#SBATCH -J RNA_Seq
#SBATCH -p batch
#SBATCH --ntasks=24
#SBATCH --mem 60gb
#SBATCH -t 8:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/RNAseq.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/RNAseq.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu


echo "This JobID for this job is ${SLURM_JOB_ID}."
echo "Done."

sacct -j $SLURM_JOB_ID --format=JobID,JobName,AllocCPUS,Elapsed,ExitCode,State,MaxRSS,TotalCPU


export LC_ALL=en_SG.utf8
export LANG=en_SG.utf8


# Author: Darrian Talamantes
# Objective: Take the output of tule annotation.smk and work with it.

# Loading Modules
module load InterProScan/5.68-100.0-foss-2022a



if [ ! -d /scratch/drt83172/Wallace_lab/RNA_SEQ/Annotation ]; then 
    mkdir -p /scratch/drt83172/Wallace_lab/RNA_SEQ/Annotation; 
fi

interproscan.sh -Xmx50g -cpu 24 -f TSV,GFF3 -goterms -b /scratch/drt83172/Wallace_lab/RNA_SEQ/Annotation/interproscan_results -i /scratch/drt83172/Wallace_lab/RNA_SEQ/transcriptome/predicted_transcripts.fasta.transdecoder_dir/predicted_transcripts.fasta.clean.transdecoder.pep -T /scratch/drt83172/Wallace_lab/RNA_SEQ/Annotation/interproscan_temp


