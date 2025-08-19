#!/bin/bash
#SBATCH -J 3b
#SBATCH -p batch	
#SBATCH --ntasks=8
#SBATCH --mem=32GB
#SBATCH -t 150:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/3b.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/3b.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu

echo "This JobID for this job is ${SLURM_JOB_ID}."
sleep 5
echo "Done."

sacct -j $SLURM_JOB_ID --format=JobID,JobName,AllocCPUS,Elapsed,ExitCode,State,MaxRSS,TotalCPU

ml Biopython/1.84-foss-2024a
module load SciPy-bundle/2024.05-gfbf-2024a

source activate GoTerm_Analysis


python 3c_gene_analysis.py