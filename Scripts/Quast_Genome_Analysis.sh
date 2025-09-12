#!/bin/bash
#SBATCH -J QUAST
#SBATCH -p batch	
#SBATCH --ntasks=8
#SBATCH --mem=32GB
#SBATCH -t 50:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/QUAST.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/QUAST.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu


ml QUAST/5.2.0

genome_dir=/scratch/drt83172/Wallace_lab/RNA_SEQ/Genome/Tall_fescue/

/scratch/drt83172/Wallace_lab/RNA_SEQ/Genome/Tall_fescue/tall_fescue_pv1.1.fasta


quast.py $genome_dir/tall_fescue_pv1.1.fasta -o $genome_dir/quast_output


