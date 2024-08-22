# Objective: This rule should go between STAR and Scallop. Once we have mapped everything to both genomes this will seperate fungal reads.

##################################################
# Getting rid of fungus
##################################################

rule grepper_big:
    input:
        big_bam = config["directories"]["star_bams"] + "Aligned.sortedByCoord.out.bam"
    conda:
        "../Conda_Envs/samtools.yaml"
    threads: 32
    output:
        filtered_bam = expand(config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.bam")
    shell:
        """
        grep -v "JAFEMN" {input.big_bam} > {output.filtered_bam}
        """

rule grepper_sep:
    input:
        sep_bams = config["directories"]["sep_bams"] + "{pairs}Aligned.sortedByCoord.out.bam"
    conda:
        "../Conda_Envs/samtools.yaml"
    threads: 8
    output:
        filtered_bams = expand(config["directories"]["filtered_bams"] + "/Aligned.sortedByCoord_filtered.out.bam")
    shell:
        """
        grep -v "JAFEMN" {input.sep_bams} > {output.filtered_bams}
        """


