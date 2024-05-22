# =================================================================================================
#     Assembling Transcriptome
# =================================================================================================

# index genome with STAR
# align reads to genome with STAR
# put bam file into trinity guided approuch 

# This resulted in an out of memory error
rule star_index:
    input:
        fasta = config["genome"],
    params:
        threads = config["params"]["star_mapping"]["threads"],
        genome_dir = config["directories"]["genome_idx"]
    conda:
        "../Conda_Envs/transcriptome.yaml"
    output:
        expand(config["directories"]["genome_idx"] + "/" + "{file}", file=star_index_files)
    shell:
        """
        STAR --runThreadN {params.threads} \
            --runMode genomeGenerate \
            --genomeDir {params.genome_dir} \
            --genomeFastaFiles {input}
        """

rule star_mapping:
    input:
        manifest = star_manifest,
        genome_dir = config["directories"]["genome_idx"]
    params:
        threads = config["params"]["star_mapping"]["threads"],
        prefix = config["directories"]["star_bams"],
    conda:
        "../Conda_Envs/transcriptome.yaml"
    output:
        bam = config["directories"]["star_bams"] + "Aligned.sortedByCoord.out.bam",
        log_out = config["directories"]["star_bams"] + "Log.out",
        log_final = config["directories"]["star_bams"] + "Log.final.out",
        sj_out = config["directories"]["star_bams"] + "SJ.out.tab"
    shell:
        """        
        STAR --runThreadN {params.threads} \
            --genomeDir {input.genome_dir} \
            --readFilesCommand zcat \
            --readFilesManifest {input.manifest} \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM SortedByCoordinate
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
