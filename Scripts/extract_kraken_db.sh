#!/bin/bash
#SBATCH -J Kdb_extract
#SBATCH -p batch
#SBATCH --ntasks=16
#SBATCH --mem 60gb
#SBATCH -t 160:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/Kdb_extract.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/Kdb_extract.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu




tar -xvzf /scratch/drt83172/Wallace_lab/RNA_SEQ/kraken_db/k2_pluspfp_20240605.tar.gz
