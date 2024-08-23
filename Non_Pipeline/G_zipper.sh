#!/bin/bash
#SBATCH -J G_zipper
#SBATCH -p batch
#SBATCH --ntasks=32
#SBATCH --mem 120gb
#SBATCH -t 140:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/G_zipper/Scripts/outfiles/RNAseq.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/G_zipper/Scripts/outfiles/RNAseq.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu

echo "This JobID for this job is ${SLURM_JOB_ID}."
sleep 5
echo "Done."

sacct -j $SLURM_JOB_ID --format=JobID,JobName,AllocCPUS,Elapsed,ExitCode,State,MaxRSS,TotalCPU




parallel/20230722-GCCcore-12.2.0

>commands.txt
while IFS= read -r file; do
    echo "gzip /scratch/drt83172/Wallace_lab/RNA_SEQ/Data/$file" >>  commands.txt
done < /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/Lists/zipme.txt
parallel --jobs 8 < commands.txt
