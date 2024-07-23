#!/bin/bash
#SBATCH -J RNA_Seq
#SBATCH -p batch
#SBATCH --ntasks=24
#SBATCH --mem 120gb
#SBATCH -t 40:00:00
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
module load eggnog-mapper/2.1.9-foss-2022a


        if [ ! -d /scratch/drt83172/Wallace_lab/RNA_SEQ/Annotation ]; then 
            mkdir -p /scratch/drt83172/Wallace_lab/RNA_SEQ/Annotation; 
        fi
emapper.py -i /scratch/drt83172/Wallace_lab/RNA_SEQ/transcriptome/predicted_transcripts.fasta.transdecoder_dir/predicted_transcripts.fasta.clean.transdecoder_subset.pep --output /scratch/drt83172/Wallace_lab/RNA_SEQ/Annotations/predicted_transcripts.emapper.annotations --cpu 24 



