# =================================================================================================
#     Assembling Transcriptome
# =================================================================================================

# index genome with STAR
# align reads to genome with STAR
# put bam file into trinity guided approuch 

rule star_index:
    input:
        fasta = config["genome"],
    params:
        threads = config["params"]["star_mapping"]["threads"],
        genome_dir = config["directories"]["genome_idx"]
    conda:
        "../Conda_Envs/transcriptome.yaml"
    threads: 32
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
        genome_files = expand(config["directories"]["genome_idx"] + "/" + "{file}", file=star_index_files)
    params:
        threads = config["params"]["star_mapping"]["threads"],
        prefix = config["directories"]["star_bams"],
        genome_dir = config["directories"]["genome_idx"]
    conda:
        "../Conda_Envs/transcriptome.yaml"
    threads: 32
    output:
        bam = config["directories"]["star_bams"] + "Aligned.sortedByCoord.out.bam",
        log_out = config["directories"]["star_bams"] + "Log.out",
        log_final = config["directories"]["star_bams"] + "Log.final.out",
        sj_out = config["directories"]["star_bams"] + "SJ.out.tab"
    shell: #ToDO: make the limitBAMsorRAM into a parameter. Right now its set to 100GB
        """        
        STAR --runThreadN {params.threads} \
            --genomeDir {params.genome_dir} \
            --readFilesCommand zcat \
            --readFilesManifest {input.manifest} \
            --outFileNamePrefix {params.prefix} \
            --limitBAMsortRAM 107089370995 \
            --outSAMtype BAM SortedByCoordinate
        """


rule star_mapping_seperate:
    input:
        fasta_fwd=trimmed + "/{pairs}R1.fq.gz",
        fasta_rev=trimmed + "/{pairs}R2.fq.gz",
        genome_files = expand(config["directories"]["genome_idx"] + "/" + "{file}", file=star_index_files)
    params:
        threads = config["params"]["star_mapping"]["threads"],
        prefix = config["directories"]["sep_bams"] + "{pairs}",
        genome_dir = config["directories"]["genome_idx"]
    conda:
        "../Conda_Envs/transcriptome.yaml"
    threads: 32
    output:
        bam = config["directories"]["sep_bams"] + "{pairs}Aligned.sortedByCoord.out.bam",
        log = config["directories"]["sep_bams"] + "{pairs}Log.out",
        log_final = config["directories"]["sep_bams"] + "{pairs}Log.final.out",
        sj = config["directories"]["sep_bams"] + "{pairs}SJ.out.tab"
    shell:
        """        
        STAR --runThreadN {params.threads} \
            --genomeDir {params.genome_dir} \
            --readFilesCommand zcat \
            --readFilesIn {input.fasta_fwd} {input.fasta_rev} \
            --outFileNamePrefix {params.prefix} \
            --limitBAMsortRAM 15000000000 \
            --outSAMtype BAM SortedByCoordinate
        """

