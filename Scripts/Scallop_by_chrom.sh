#!/bin/bash
#SBATCH -J Scallop_Split
#SBATCH -p highmem_p	
#SBATCH --ntasks=32
#SBATCH --mem=500GB
#SBATCH -t 160:00:00
#SBATCH --output=/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/Scallop_Split.%j.out
#SBATCH -e /scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/outfiles/Scallop_Split.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user drt83172@uga.edu

echo "This JobID for this job is ${SLURM_JOB_ID}."
sleep 5
echo "Done."

sacct -j $SLURM_JOB_ID --format=JobID,JobName,AllocCPUS,Elapsed,ExitCode,State,MaxRSS,TotalCPU






# Module loading
module load Anaconda3/2022.10
module load snakemake/7.22.0-foss-2022a
module load parallel/20240322-GCCcore-13.2.0
source activate smatools


#!/bin/bash

set -e  # Exit on error
set -o pipefail  # Catch pipeline errors


# Define directories (update these based on your config)
GENOME="/scratch/drt83172/Wallace_lab/RNA_SEQ/Genome/Tall_fescue/tall_fescue_pv1.1.fasta"
FILTERED_BAM_DIR="/scratch/drt83172/Wallace_lab/RNA_SEQ/filtered_bams/big"
BIG_BAM_CHROM="/scratch/drt83172/Wallace_lab/RNA_SEQ/filtered_bams/big/split"
SCALLOP_OUT="/scratch/drt83172/Wallace_lab/RNA_SEQ/transcriptome"
FINAL_GTF="/scratch/drt83172/Wallace_lab/RNA_SEQ/transcriptome_final/Fescue_transcriptome.gtf"

# 1️⃣ **Fix BAM Header**
echo "Fixing BAM header..."
cd "$FILTERED_BAM_DIR"
samtools faidx "$GENOME"
cut -f1,2 "${GENOME}.fai" | awk '{print "@SQ\tSN:"$1"\tLN:"$2}' > new_header.sam
samtools reheader new_header.sam Aligned.sortedByCoord_filtered.out.bam  > Aligned.sortedByCoord_filtered_fixed.out.bam

# 2️⃣ **Index the Fixed BAM**
echo "Indexing BAM..."
samtools index -c Aligned.sortedByCoord_filtered_fixed.out.bam

# 3️⃣ **Split BAM by Chromosome**
echo "Splitting BAM by chromosome..."
mkdir -p "$BIG_BAM_CHROM"
samtools idxstats Aligned.sortedByCoord_filtered_fixed.out.bam | cut -f1 | grep -v '*' | while read -r chr; do
    echo "Processing chromosome: $chr"
    samtools view -b -h Aligned.sortedByCoord_filtered_fixed.out.bam "$chr" | samtools sort -o "$BIG_BAM_CHROM/${chr}.bam"
    samtools index -c "$BIG_BAM_CHROM/${chr}.bam"
done

# 4️⃣ **Run Scallop2 in Parallel**
echo "Running Scallop2 in parallel..."
mkdir -p "$SCALLOP_OUT"

export SCALLOP_OUT  # Ensure the variable is available for parallel

parallel --jobs 8 --halt now,fail=1 '
    source activate scallop2 &&
    chrom=$(basename {} .bam) &&
    echo "Assembling transcripts for $chrom..." &&
    scallop2 --num-threads 8 -i {} -o "$SCALLOP_OUT/${chrom}.gtf" &&
    conda deactivate
' ::: "$BIG_BAM_CHROM"/*.bam

# 5️⃣ **Merge GTF Files**
echo "Merging GTF files..."
cat "$SCALLOP_OUT"/*.gtf | grep -v '^#' | sort -k1,1 -k4,4n > "$FINAL_GTF"

echo "Pipeline completed successfully!"







# Aligned.sortedByCoord_filtered
# test_out.bam
