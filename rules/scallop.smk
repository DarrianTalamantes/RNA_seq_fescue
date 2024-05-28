rule scallop2:
    input:
        bam= config["directories"]["star_bams"] + "Aligned.sortedByCoord.mapped.out.bam"
    conda:
        "../Conda_Envs/scallop2.yaml"
    output:
        gtf=config["scallop"]["output_file"]
    shell:
        """
        scallop2 -i {input.bam} -o {output.gtf}
        """


