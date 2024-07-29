#!/bin/bash
#SBATCH -J RNA_Seq
#SBATCH -p batch
#SBATCH --ntasks=12
#SBATCH --mem 60gb
#SBATCH -t 5:00:00
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
module load BBMap/39.01-GCC-12.2.0


stats.sh /scratch/drt83172/Wallace_lab/RNA_SEQ/Genome/tall_fescuev0.1.fa
