#!/bin/bash
#SBATCH -J Fixer
#SBATCH -p batch
#SBATCH --ntasks=4
#SBATCH --mem 100gb
#SBATCH -t 140:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/Fixer.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/Fixer.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu

echo "This JobID for this job is ${SLURM_JOB_ID}."
sleep 5
echo "Done."

ml STAR/2.7.10b-GCC-11.3.0     
    
        STAR --runThreadN 4 \
            --genomeDir /scratch/drt83172/Wallace_lab/RNA_SEQ/Genome_idx \
            --readFilesCommand zcat \
            --readFilesIn /scratch/drt83172/Wallace_lab/RNA_SEQ/Trimmed/CTE46P_O17_HP3_S356_F59_R1.fq.gz /scratch/drt83172/Wallace_lab/RNA_SEQ/Trimmed/CTE46P_O17_HP3_S356_F59_R2.fq.gz \
            --outFileNamePrefix /scratch/drt83172/Wallace_lab/RNA_SEQ/sep_bams/CTE46P_O17_HP3_S356_F59_ \
            --limitBAMsortRAM 15000000000 \
            --outSAMtype BAM SortedByCoordinate
