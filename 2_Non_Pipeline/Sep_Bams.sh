#!/bin/bash
#SBATCH -J sepStar
#SBATCH -p batch
#SBATCH --ntasks=32
#SBATCH --mem 120gb
#SBATCH -t 150:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/sepStar.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/sepStar.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu

echo "This JobID for this job is ${SLURM_JOB_ID}."
sleep 5
echo "Done."

sacct -j $SLURM_JOB_ID --format=JobID,JobName,AllocCPUS,Elapsed,ExitCode,State,MaxRSS,TotalCPU

module load Anaconda3/2022.10
module load snakemake/6.9.1-Mamba-4.11.0-4

# Line allows you to use conda
. $(conda info --root)/etc/profile.d/conda.sh 

source activate transcriptome

export LC_ALL=en_SG.utf8
export LANG=en_SG.utf8

for sample in $(cat /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/paired_list.txt)
    do
        STAR --runThreadN 24 \
            --genomeDir /scratch/drt83172/Wallace_lab/RNA_SEQ/Genome_idx \
            --readFilesIn  /scratch/drt83172/Wallace_lab/RNA_SEQ/kraken/non_fungal/${sample}R1.fq /scratch/drt83172/Wallace_lab/RNA_SEQ/kraken/non_fungal/${sample}R2.fq \
            --outFileNamePrefix /scratch/drt83172/Wallace_lab/RNA_SEQ/sep_bams/$sample \
            --limitBAMsortRAM 15000000000 \
            --outSAMtype BAM SortedByCoordinate
    done