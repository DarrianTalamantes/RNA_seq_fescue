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
        chunk_prefix= config["directories"]["filtered_bam_big"] + "/chunk_",
        filtered_sam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.sam",
        inter_sam = config["directories"]["big_bam"] + "Aligned.sortedByCoord.out.sam",
        lines_per_chunk = config["fungal_removal"]["chunk"]
    threads: 32
    output:
        filtered_bam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.bam"
    shell:
        """
        # Create output directory if it doesn't exist
        if [ ! -d {params.output_dir} ]; then
            mkdir -p {params.output_dir}
        fi

        # Convert BAM to SAM
        samtools view {input.big_bam} > {params.inter_sam}

        # Split the SAM file into chunks
        split -l {params.lines_per_chunk} {params.inter_sam} {params.chunk_prefix}

        # Use parallel to grep and filter each chunk
        ls {params.chunk_prefix}* | parallel -j {threads} "grep -v 'JAFEMN' {{}} > {{}}.out"

        # Concatenate filtered chunks into a single SAM file
        cat {params.chunk_prefix}*.out > {params.filtered_sam}

        # Clean up intermediate files
        rm {params.chunk_prefix}*
        
        # Convert filtered SAM back to BAM
        samtools view -b {params.filtered_sam} > {output.filtered_bam}

        # Remove intermediate SAM files
        rm {params.filtered_sam}
        rm {params.inter_sam}
        """

        
rule grepper_sep:
    input:
        sep_bams = config["directories"]["sep_bams"] + "{pairs}Aligned.sortedByCoord.out.bam"
    conda:
        "../Conda_Envs/samtools.yaml"
    params:
        output_dir = config["directories"]["filtered_bams"]
    threads: 8
    output:
        filtered_bams = config["directories"]["filtered_bams"] + "/{pairs}Aligned.sortedByCoord_filtered.out.bam"
    shell:
        """
        # Create output directory if it doesn't exist
        if [ ! -d {params.output_dir} ]; then
            mkdir -p {params.output_dir}
        fi

        # Convert BAM to SAM, filter with grep, and save the filtered SAM
        samtools view -@ {threads} {input.sep_bams} | grep -v "JAFEMN" | samtools view -b -@ {threads} - > {output.filtered_bams}
        """



