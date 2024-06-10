        STAR --runThreadN 4 \
            --genomeDir {/scratch/drt83172/Wallace_lab/RNA_SEQ/Genome_idx} \
            --readFilesCommand zcat \
            --readFilesIn /scratch/drt83172/Wallace_lab/RNA_SEQ/Trimmed/CTE46P_O17_HP3_S356_F59_R1.fq.gz /scratch/drt83172/Wallace_lab/RNA_SEQ/Trimmed/CTE46P_O17_HP3_S356_F59_R2.fq.gz \
            --outFileNamePrefix /scratch/drt83172/Wallace_lab/RNA_SEQ/sep_bams/CTE46P_O17_HP3_S356_F59_ \
            --limitBAMsortRAM 15000000000 \
            --outSAMtype BAM SortedByCoordinate
