#!/bin/bash
#SBATCH -J FastQC
#SBATCH -p batch
#SBATCH --ntasks=12
#SBATCH --mem 40gb
#SBATCH -t 140:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/FastQC.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/FastQC.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu

echo "This JobID for this job is ${SLURM_JOB_ID}."
sleep 5
echo "Done."

sacct -j $SLURM_JOB_ID --format=JobID,JobName,AllocCPUS,Elapsed,ExitCode,State,MaxRSS,TotalCPU

module load FastQC/0.11.9-Java-11

./scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/Scripts/RunFastQC.sh /scratch/drt83172/Wallace_lab/RNA_SEQ/Trimmed
