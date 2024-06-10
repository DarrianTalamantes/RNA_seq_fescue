        STAR --runThreadN {4} \
            --genomeDir {/scratch/drt83172/Wallace_lab/RNA_SEQ/Genome/tall_fescuev0.1.fa} \
            --readFilesCommand zcat \
            --readFilesIn CTE46P_O17_HP3_S356_F59_R1.fq.gz CTE46P_O17_HP3_S356_F59_R2.fq.gz \
            --outFileNamePrefix CTE46P_O17_HP3_S356_F59_ \
            --limitBAMsortRAM 15000000000 \
            --outSAMtype BAM SortedByCoordinate
