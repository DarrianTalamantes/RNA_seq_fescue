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

# rule star_pe_multi:
#     input:
#         # We actually need manifest file here not whatever these inputs are
#         manifest = star_manifest,
#         genome_index=config["directories"]["genome_idx"],
#     output:
#         # see STAR manual for additional output files
#         aln= align + "{pairs}/pe_aligned.sam",
#         log= align + "{pairs}/Log.out",
#         sj= align + "{pairs}/SJ.out.tab",
#         unmapped=[ align + "{pairs}/unmapped.1.fastq.gz", align + "{pairs}/unmapped.2.fastq.gz"],
#     conda:
#         "Conda_Envs/transcriptome.yaml"
#     threads: 8
#     params:
#         threads = 8
#         prefix = align
#     shell:
#         """
#     STAR --genomeDir {input.genome_index} \
#     --readFilesManifest {input.manifest} \
#     --runThreadN {params.threads} \
#     --outSAMtype BAM SortedByCoordinate \
#     --outFileNamePrefix {params.prefix}/aligned_reads
#         """




# rule star_alignment:
#     input:
#         reads1 = trimmed + "/{pairs}R1.fq.gz",
#         reads2 = trimmed + "/{pairs}R2.fq.gz",
#         genome_index = config["directories"]["genome_idx"]
#     output:
#         alignment = align + "aligned_reads.bam"
#     conda:
#         'Conda_Envs/transcriptome.yaml'
#     threads: 10
#     shell:
#     """
#     STAR --genomeDir {input.genome_index} \
#     --readFilesIn {input.reads1} {input.reads2} \
#     --runThreadN {params.threads} \
#     --outSAMtype BAM SortedByCoordinate \
#     --outFileNamePrefix aligned_reads
#     """


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
