# Objective: This rule should go between STAR and Scallop. Once we have mapped everything to both genomes this will seperate fungal reads.

##################################################
# Getting rid of fungus
##################################################

rule filter_epichloe_out:
    input:
        big_bam=config["directories"]["big_bam"] + "Aligned.sortedByCoord.out.bam"
    output:
        filtered_bam=config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.bam"
    threads: config["fungal_removal"]["threads"]
    shell:
        """
        # Extract @RG lines and keep only those related to JAFEMN
        samtools view -H {input.big_bam} | grep "^@RG" | grep "JAFEMN" | cut -f 2 | sed 's/ID://g' > keep_rg.txt

        # Filter the BAM file to include only those read groups
        samtools view --threads {threads} -b -h $(cat keep_rg.txt | xargs -I {{}} echo -r {{}}) {input.big_bam} > {output.filtered_bam}
        """

rule grepper_sep:
    input:
        sep_bams = config["directories"]["sep_bams"] + "{pairs}Aligned.sortedByCoord.out.bam"
    conda:
        "../../Conda_Envs/samtools.yaml"
    params:
        output_dir = config["directories"]["filtered_bams"]
    threads: 8
    output:
        filtered_bams = config["directories"]["filtered_bams"] + "/{pairs}Aligned.sortedByCoord_filtered.out.bam"
    log:
        temp("logs/{pairs}_grepper.log")
    shell:
        """
        # Create output directory if it doesn't exist
        if [ ! -d {params.output_dir} ]; then
            mkdir -p {params.output_dir}
        fi

        # Step 1: Convert BAM to SAM and check if it's successful
        samtools view -h -@ {threads} {input.sep_bams} > {params.output_dir}/{wildcards.pairs}.tmp.sam 2>> {log}
        if [ $? -ne 0 ]; then
            echo "Error in converting BAM to SAM" >> {log}
            exit 1
        fi

        # Step 2: Preserve the header and filter the SAM using grep (only for alignments)
        grep "^@" {params.output_dir}/{wildcards.pairs}.tmp.sam > {params.output_dir}/{wildcards.pairs}.header.sam
        grep "JAFEMN" {params.output_dir}/{wildcards.pairs}.tmp.sam | grep -v "^@" >> {params.output_dir}/{wildcards.pairs}.header.sam
        mv {params.output_dir}/{wildcards.pairs}.header.sam {params.output_dir}/{wildcards.pairs}.filtered.sam

        if [ $? -ne 0 ] || [ ! -s {params.output_dir}/{wildcards.pairs}.filtered.sam ]; then
            echo "Error in filtering SAM or no output generated" >> {log}
            exit 1
        fi

        # Step 3: Convert filtered SAM back to BAM
        samtools view -b -@ {threads} {params.output_dir}/{wildcards.pairs}.filtered.sam > {output.filtered_bams} 2>> {log}
        if [ $? -ne 0 ]; then
            echo "Error in converting filtered SAM to BAM" >> {log}
            exit 1
        fi

        # Clean up temporary files
        rm {params.output_dir}/{wildcards.pairs}.tmp.sam {params.output_dir}/{wildcards.pairs}.filtered.sam
        """





