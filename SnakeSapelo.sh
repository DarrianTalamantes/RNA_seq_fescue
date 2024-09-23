#!/bin/bash
#SBATCH -J FungalBig
#SBATCH -p batch
#SBATCH --ntasks=32
#SBATCH --mem 120gb
#SBATCH -t 150:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/FungalBig.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/FungalBig.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu

echo "This JobID for this job is ${SLURM_JOB_ID}."
sleep 5
echo "Done."

sacct -j $SLURM_JOB_ID --format=JobID,JobName,AllocCPUS,Elapsed,ExitCode,State,MaxRSS,TotalCPU

module load Anaconda3/2022.10
module load snakemake/6.9.1-Mamba-4.11.0-4
source activate snakemake

export LC_ALL=en_SG.utf8
export LANG=en_SG.utf8

snakemake --use-conda --cores 32  -s RNAseq.smk --verbose --rerun-triggers mtime

# --rerun-incomplete