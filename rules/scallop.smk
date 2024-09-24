rule scallop2:
    input:
        bam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.bam"
    conda:
        "../Conda_Envs/scallop2.yaml"
    output:
        gtf = config["scallop"]["output_file"]
    shell:
        """
        scallop2 -i {input.bam} -o {output.gtf}
        """


