#!/bin/bash
#SBATCH -J EnTAP
#SBATCH -p batch
#SBATCH --ntasks=4
#SBATCH --mem 60gb
#SBATCH -t 30:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/RNA_SEQ/EnTAP/slurm_outputs/EnTAP.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/RNA_SEQ/EnTAP/slurm_outputs/EnTAP.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu


sacct -j $SLURM_JOB_ID --format=JobID,JobName,AllocCPUS,Elapsed,ExitCode,State,MaxRSS,TotalCPU

# Loading Modules
ml EnTAP/2.2.0-GCCcore-11.3.0

# Code to be Ran

EnTAP --run --run-ini /scratch/drt83172/Wallace_lab/RNA_SEQ/EnTAP/Params/entap_run.params --entap-ini /scratch/drt83172/Wallace_lab/RNA_SEQ/EnTAP/Params/entap_config_ini

