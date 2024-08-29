# =================================================================================================
#     Assembling Transcriptome
# =================================================================================================

# This should go after kraken
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







# This uses a manifest file to list everything that should be mapped into the one big file
#Note: This manifest file has to change to the new kraken stuff.
if config["use_ignored_rule"]:
    rule star_mapping:
        input:
            genome_files = expand(config["directories"]["genome_idx"] + "/" + "{file}", file=star_index_files)
        params:
            threads = config["star_mapping"]["threads"],
            prefix = config["directories"]["star_bams"],
            genome_dir = config["directories"]["genome_idx"],
            manifest = config["star_mapping"]["star_manifest"]

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
            ls        
            STAR --runThreadN {params.threads} \
                --genomeDir {params.genome_dir} \
                --readFilesCommand zcat \
                --readFilesManifest {params.manifest} \
                --outFileNamePrefix {params.prefix} \
                --limitBAMsortRAM 107089370995 \
                --outSAMtype BAM SortedByCoordinate
            """











# This rule maps all the files seperatly with many different output files.
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
######## This is a a rough draft of the new version of this rule
rule star_mapping:
    input:
        krakened_fq = config["kraken"]["classified"] + "/krakened_{pairs}.fq.gz",
        star_index = config["star"]["index"]  # Path to the STAR index
    output:
        aligned_bam = config["star"]["output"] + "/{pairs}.Aligned.sortedByCoord.out.bam"
    params:
        threads = config["star"]["threads"],
        out_prefix = config["star"]["output"] + "/{pairs}."  # Prefix for STAR output files
    conda:
        "../Conda_Envs/star.yaml"
    shell:
        """
        STAR --runThreadN {params.threads} \
             --genomeDir {input.star_index} \
             --readFilesIn {input.krakened_fq} \
             --readFilesCommand zcat \
             --outFileNamePrefix {params.out_prefix} \
             --outSAMtype BAM SortedByCoordinate
        """