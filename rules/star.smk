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

rule seperate_mapped:
    input:
        config["directories"]["star_bams"] + "Aligned.sortedByCoord.out.bam"
    conda:
        "../Conda_Envs/samtools.yaml"
    output:
        config["directories"]["star_bams"] + "Aligned.sortedByCoord.mapped.out.bam"
    shell:
        """
        samtools view -b -F 4 {input} > {output}
        """


rule bam_seperation:
    input:
        bam = config["directories"]["star_bams"] + "Aligned.sortedByCoord.mapped.out.bam"
    conda:
        "../Conda_Envs/samtools.yaml"
    output:
        expand(sep_bams + "/{sample}.bam", sample=SAMPLES) 
    shell:
        """
            for sample in {SAMPLES}; do
                samtools view -b -r $sample {input.bam} > {sep_bams}/$sample.bam;
            done
        """

