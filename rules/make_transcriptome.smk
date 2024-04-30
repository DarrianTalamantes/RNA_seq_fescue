# =================================================================================================
#     Assembling Transcriptome
# =================================================================================================

# index genome with STAR
# align reads to genome with STAR
# put bam file into trinity guided approuch 

rule star_index:
    input:
        genome_directory = genomedic
        genome_file =trimmed + genome
    output:
        directory(genomedic)
    conda:
        'Conda_Envs/transcriptome.yaml',
    threads: 10
    shell:
    """
    STAR --runThreadN {threads} \
    --runMode genomeGenerate \
    --genomeDir {input.genome_directory} \
    --genomeFastaFiles {input.genome_file} \
    --sjdbOverhang 101
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

'''
