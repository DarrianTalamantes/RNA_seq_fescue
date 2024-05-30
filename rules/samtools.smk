# Objective: These rules will seperate the bam file into many smaller files.


#############################
# This is deactivated cause I think it is an extra step that might not do anything
#############################
# rule seperate_mapped:
#     input:
#         config["directories"]["star_bams"] + "Aligned.sortedByCoord.out.bam"
#     conda:
#         "../Conda_Envs/samtools.yaml"
#     output:
#         config["directories"]["star_bams"] + "Aligned.sortedByCoord.mapped.out.bam"
#     shell:
#         """
#         samtools view -b -F 4 {input} > {output}
#         """


rule bam_seperation:
    input:
        bam = config["directories"]["star_bams"] + "Aligned.sortedByCoord.out.bam"
    conda:
        "../Conda_Envs/samtools.yaml"
    output:
        expand(sep_bams + "{pairs}.bam", pairs=PAIRS) 
    params:
        dir = config["directories"]["sep_bams"],
        logsheet = config["params"]["star_mapping"]["log"],
        log_dir = config["directories"]["log"]
    shell:
        """
        if [ ! -d {params.log_dir} ]; then \
            mkdir -p {params.log_dir}; fi

        for pair in {PAIRS}; do
            echo "Processing $pair" >> {params.logsheet}
            samtools view -b -r ${{pair}}R1 {input.bam} > {params.dir}${{pair}}.bam
        done
        """