#!/bin/bash
#SBATCH -J EnTAP
#SBATCH -p batch
#SBATCH --ntasks=16
#SBATCH --mem 80gb
#SBATCH -t 160:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/RNA_SEQ/EnTAP/slurm_outputs/EnTAP.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/RNA_SEQ/EnTAP/slurm_outputs/EnTAP.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu


sacct -j $SLURM_JOB_ID --format=JobID,JobName,AllocCPUS,Elapsed,ExitCode,State,MaxRSS,TotalCPU

# Loading Modules
ml EnTAP/2.2.0-GCCcore-11.3.0

# # Must run this first with refseq_plants in the .params --database, then put path to dimond database in the --database .params file
# EnTAP --config --run-ini /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/EnTAP_scripts/entap_run.params --entap-ini /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/EnTAP_scripts/entap_config.ini

# Code to be Ran
EnTAP --run --run-ini /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/EnTAP_scripts/entap_run.params --entap-ini /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/EnTAP_scripts/entap_config.ini

