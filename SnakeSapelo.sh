#!/bin/bash
#SBATCH -J Scallop_Split
#SBATCH -p highmem	
#SBATCH --ntasks=32
#SBATCH --mem=500GB
#SBATCH -t 160:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/Scallop_Split.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/Scallop_Split.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu

echo "This JobID for this job is ${SLURM_JOB_ID}."
sleep 5
echo "Done."

sacct -j $SLURM_JOB_ID --format=JobID,JobName,AllocCPUS,Elapsed,ExitCode,State,MaxRSS,TotalCPU

module load Anaconda3/2022.10
module load snakemake/7.22.0-foss-2022a
source activate snakemake

export LC_ALL=en_SG.utf8
export LANG=en_SG.utf8

snakemake --use-conda --cores 32 -s RNAseq.smk --verbose --rerun-incomplete --forcerun fix_bam_header

#--rerun-triggers mtime 
# --rerun-incomplete