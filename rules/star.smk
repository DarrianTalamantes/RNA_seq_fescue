# =================================================================================================
#     Assembling Transcriptome
# =================================================================================================

# This should go after kraken
# index genome with STAR
# align reads to genome with STAR
# In this usaage I have combined the epi and Fescue genome and align the reads to both (will be seperated later) 

rule star_index:
    input:
        fasta = config["combined_genome"],
    params:
        threads = config["star_mapping"]["threads"],
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

rule create_star_manifest:
    input:
        # Using the outputs from the Kraken rule
        fwd_files = expand(config["kraken"]["non_fungal"] + "/{pairs}R1.fq", pairs=PAIRS),
        rev_files = expand(config["kraken"]["non_fungal"] + "/{pairs}R2.fq", pairs=PAIRS)
    output:
        manifest = config["star_mapping"]["star_manifest"]
    run:
        import os

        # Create a list to hold lines of the manifest
        manifest_lines = []

        # Iterate over each pair of forward and reverse FASTQ files
        for fwd, rev in zip(input.fwd_files, input.rev_files):
            sample_name = os.path.basename(fwd).replace("_R1.fq", "")
            manifest_lines.append(f"{fwd}\t{rev}\t{sample_name}\n")

        # Write the lines to the manifest file
        with open(output.manifest, "w") as manifest_file:
            manifest_file.writelines(manifest_lines)





# This uses a manifest file to list everything that should be mapped into the one big file
#Note: This manifest file has to change to the new kraken stuff.
if config["use_ignored_rule"]:
    rule star_mapping:
        input:
            genome_files = expand(config["directories"]["genome_idx"] + "/" + "{file}", file=star_index_files),
            manifest = config["star_mapping"]["star_manifest"]
        params:
            threads = config["star_mapping"]["threads"],
            prefix = config["directories"]["big_bam"],
            genome_dir = config["directories"]["genome_idx"]

        conda:
            "../Conda_Envs/transcriptome.yaml"
        threads: config["star_mapping"]["threads"]
        output:
            bam = config["directories"]["big_bam"] + "Aligned.sortedByCoord.out.bam",
            log_out = config["directories"]["big_bam"] + "Log.out",
            log_final = config["directories"]["big_bam"] + "Log.final.out",
            sj_out = config["directories"]["big_bam"] + "SJ.out.tab"
        shell: #ToDO: make the limitBAMsorRAM into a parameter. Right now its set to 300GB
            """ 
        if [ ! -d {params.prefix} ]; then
            mkdir -p {params.prefix}; 
        fi
            STAR --runThreadN {params.threads} \
                --genomeDir {params.genome_dir} \
                --readFilesManifest {input.manifest} \
                --outFileNamePrefix {params.prefix} \
                --limitBAMsortRAM 307089370995 \
                --outSAMtype BAM SortedByCoordinate
            """

# This rule maps all the files seperatly with many different output files.
# rule star_mapping_seperate:
    input:
        extracted_fwd = config["kraken"]["non_fungal"] + "/{pairs}R1.fq",
        extracted_rev = config["kraken"]["non_fungal"] + "/{pairs}R2.fq",
        genome_files = expand(config["directories"]["genome_idx"] + "/" + "{file}", file=star_index_files)
    params:
        threads = config["star_mapping"]["threads_sep"],
        prefix = config["directories"]["sep_bams"] + "{pairs}",
        genome_dir = config["directories"]["genome_idx"]
    conda:
        "../Conda_Envs/transcriptome.yaml"
    threads: config["star_mapping"]["threads_sep"]
    output:
        bam = config["directories"]["sep_bams"] + "{pairs}Aligned.sortedByCoord.out.bam",
        log = config["directories"]["sep_bams"] + "{pairs}Log.out",
        log_final = config["directories"]["sep_bams"] + "{pairs}Log.final.out",
        sj = config["directories"]["sep_bams"] + "{pairs}SJ.out.tab"
    shell:
        """        
        STAR --runThreadN {params.threads} \
            --genomeDir {params.genome_dir} \
            --readFilesIn {input.extracted_fwd} {input.extracted_rev} \
            --outFileNamePrefix {params.prefix} \
            --limitBAMsortRAM 15000000000 \
            --outSAMtype BAM SortedByCoordinate
        """
