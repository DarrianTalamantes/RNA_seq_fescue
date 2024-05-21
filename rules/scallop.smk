rule scallop2:
    input:
        bam= config["directories"]["star_bams"] + "Aligned.sortedByCoord.out.bam"
    params:
        path="3_STAR/3b_map"
    conda:
        "../Conda_Envs/transcriptome.yaml"
    output:
        "{sample}"
    shell:
        """
        scallop2 -i {params.path}/{input.bam} -o {output}.gtf
        """


