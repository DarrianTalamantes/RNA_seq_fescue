# =================================================================================================
#     Assembling Transcriptome
# =================================================================================================

# index genome with STAR
# align reads to genome with STAR
# put bam file into trinity guided approuch 

# This resulted in an out of memory error
rule star_index:
    input:
        fasta = genome,
    output:
        directory(config["directories"]["genome_idx"])
    threads:10
    params:
        extra = ""
    wrapper:
        "0.49.0/bio/star/index"

rule star_mapping:
    input:
        manifest = star_manifest,
        genome_dir = config["directories"]["genome_idx"]
    params:
        threads = config["params"]["star_mapping"]["threads"],
        compcomm = config["params"]["star_mapping"]["compcomm"]
    output:
        bam = config["directories"]["star_bams"] + "Aligned.sortedByCoord.out.bam"
    run:
        shell("STAR --runThreadN {params.threads} \
            --genomeDir {input.genome_dir} \
            --readFilesCommand zcat \
            --readFilesManifest {input.manifest} \
            --outFileNamePrefix {output.bam} \
            --outSAMtype BAM SortedByCoordinate")





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
