# Objective: This rule should go between STAR and Scallop. Once we have mapped everything to both genomes this will seperate fungal reads.

##################################################
# Getting rid of fungus
##################################################

rule grepper_big:
    input:
        big_bam = config["directories"]["big_bam"] + "Aligned.sortedByCoord.out.bam"
    conda:
        "../Conda_Envs/samtools.yaml"
    params:
        output_dir= config["directories"]["filtered_bam_big"],
        chunk= config["fungal_removal"]["chunk"],
        filtered_sam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.sam"
    threads: 32
    output:
        filtered_bam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.bam"
    shell:
        """
        if [ ! -d {params.output_dir} ]; then
            mkdir -p {params.output_dir}; 
        fi

        samtools view {input.big_bam} | split -l 10000000 chunk_
        ls chunk_* | parallel -j 8 "grep -v 'JAFEMN' {{}} > {{}}.out"
        cat chunk_*.out > {params.filtered_sam}
        rm chunk_*
        samtools view -b {params.filtered_sam} > {output.filtered_bam}
        rm {params.filtered_sam}
        """

rule grepper_sep:
    input:
        sep_bams = config["directories"]["sep_bams"] + "{pairs}Aligned.sortedByCoord.out.bam"
    conda:
        "../Conda_Envs/samtools.yaml"
    params:
        output_dir = config["directories"]["filtered_bams"],
        filtered_sam = config["directories"]["filtered_bams"] + "/{pairs}Aligned.sortedByCoord_filtered.out.sam"
    threads: 8
    output:
        filtered_bams = config["directories"]["filtered_bams"] + "/{pairs}Aligned.sortedByCoord_filtered.out.bam"
    shell:
        """
        if [ ! -d {params.output_dir} ]; then
            mkdir -p {params.output_dir}; 
        fi

        samtools view {input.sep_bams} | grep -v "JAFEMN" > {params.filtered_sam}
        samtools view -b {params.filtered_sam} > {output.filtered_bams}
        rm {params.filtered_sam} 
        """


