rule scallop2:
    input:
        bam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.bam"
    conda:
        "../Conda_Envs/scallop2.yaml"
    params:
        threads = config["scallop"]["threads"]
    threads: config["scallop"]["threads"]
    output:
        gtf = config["scallop"]["output_file"]
    shell:
        """
        scallop2 --num-threads {params.threads} -i {input.bam} -o {output.gtf}
        """


