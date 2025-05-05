#!/bin/bash
#SBATCH -J gffread
#SBATCH -p batch
#SBATCH --ntasks=8
#SBATCH --mem 32gb
#SBATCH -t 30:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/RNA_SEQ/EnTAP/slurm_outputs/gffread.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/RNA_SEQ/EnTAP/slurm_outputs/gffread.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu


sacct -j $SLURM_JOB_ID --format=JobID,JobName,AllocCPUS,Elapsed,ExitCode,State,MaxRSS,TotalCPU

# Loading Modules
ml gffread/0.12.7-GCCcore-11.3.0

# Code to be Ran

# This script will take the Scallop output and create a .fa file for use with EnTaP

scallop_output="/scratch/drt83172/Wallace_lab/RNA_SEQ/transcriptome/big/Fescue_transcriptome.gtf"
genome="/scratch/drt83172/Wallace_lab/RNA_SEQ/Genome/Tall_fescue/tall_fescue_pv1.1.fasta"
transcripts="/scratch/drt83172/Wallace_lab/RNA_SEQ/EnTAP/Transcripts/fescue_transcripts.fasta"
gffread $scallop_output -g $genome -w $transcripts -F 
