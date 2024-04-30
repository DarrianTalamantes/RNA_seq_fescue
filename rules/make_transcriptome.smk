# =================================================================================================
#     Assembling Transcriptome
# =================================================================================================

# index genome with STAR
# align reads to genome with STAR
# put bam file into trinity guided approuch 

rule star_index:
    input:
        genome_file =trimmed + genome
    output:
        directory(config["directories"]["genome_idx"])
    conda:
        'Conda_Envs/transcriptome.yaml',
    threads: 10
    shell:
    """
    STAR --runThreadN {threads} \
    --runMode genomeGenerate \
    --genomeDir {output} \
    --genomeFastaFiles {input.genome_file} \
    --sjdbOverhang 101
    """

rule star_alignment:
    input:
        reads1 = trimmed + "/{pairs}R1.fq.gz",
        reads2 = trimmed + "/{pairs}R2.fq.gz",
        genome_index=directory("genome_index")
    output:
        alignment = align + "aligned_reads.bam"
    conda:
        'Conda_Envs/transcriptome.yaml'
    threads: 10
    shell:
    """
    STAR --genomeDir {input.genome_index} \
    --readFilesIn {input.reads1} {input.reads2} \
    --runThreadN {params.threads} \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix aligned_reads

    """


'''
Trinity --seqType fq --left reads_1.fq --right reads_2.fq --CPU 6 --max_memory 20G 

This is only for 1 sample at a time
STAR --runThreadN [num threads] \
        --genomeDir [location of indexed genome] \
        --readFilesCommand gunzip -c \
        --readFilesIn [forward read] [reverse read] \
        --outFileNamePrefix [prefix for bam files and other files] \
        --outSAMtype BAM SortedByCoordinate





rule align_reads_with_star:
    input:
        reads1="reads_R1.fastq",
        reads2="reads_R2.fastq",
        genome_index=directory("genome_index")  # Use the output of the indexing rule as input
    output:
        "aligned_reads.bam"
    params:
        threads=8
    shell:
        """
        STAR --genomeDir {input.genome_index} \
                --readFilesIn {input.reads1} {input.reads2} \
                --runThreadN {params.threads} \
                --outSAMtype BAM SortedByCoordinate \
                --outFileNamePrefix aligned_reads
        """

'''
