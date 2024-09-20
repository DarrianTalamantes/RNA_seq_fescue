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
        output_dir = config["directories"]["filtered_bam_big"],
        chunk_prefix = config["directories"]["filtered_bam_big"] + "/chunk_",
        filtered_sam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.sam",
        inter_sam = config["directories"]["big_bam"] + "Aligned.sortedByCoord.out.sam",
        lines_per_chunk = config["fungal_removal"]["lines_per_chunk"]
    threads: 32
    output:
        filtered_bam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.bam"
    log:
        temp("logs/grepper_big.log")
    shell:
        """
        echo "Starting grepper_big rule" >> {log}
        
        # Create output directory if it doesn't exist
        if [ ! -d {params.output_dir} ]; then
            mkdir -p {params.output_dir}
            echo "Created output directory: {params.output_dir}" >> {log}
        fi

        # Step 1: Convert BAM to SAM with header (-h)
        samtools view -h {input.big_bam} > {params.inter_sam} 2>> {log}
        if [ $? -ne 0 ]; then
            echo "Error in converting BAM to SAM" >> {log}
            exit 1
        fi
        echo "Converted BAM to SAM successfully" >> {log}

        # Step 2: Extract the SAM header
        grep "^@" {params.inter_sam} > {params.output_dir}/sam_header.sam
        if [ $? -ne 0 ]; then
            echo "Error extracting SAM header" >> {log}
            exit 1
        fi
        echo "Extracted SAM header successfully" >> {log}

        # Step 3: Split the SAM file (without the header) into chunks
        grep -v "^@" {params.inter_sam} | split -l {params.lines_per_chunk} - {params.chunk_prefix}
        if [ $? -ne 0 ]; then
            echo "Error splitting SAM file" >> {log}
            exit 1
        fi
        echo "Split SAM file into chunks successfully" >> {log}

        # Check if any chunks were created
        if [ -z "$(ls {params.chunk_prefix}*)" ]; then
            echo "No chunks were created during splitting" >> {log}
            exit 1
        fi

        # Step 4: Use parallel to grep and filter each chunk
        echo "Filtering chunks..." >> {log}
        ls {params.chunk_prefix}* | parallel -j {threads} "grep -v 'JAFEMN' {{}} > {{}}.out"

        # Check if the filtered chunks are generated
        if [ -z "$(ls {params.chunk_prefix}*.out)" ]; then
            echo "No filtered chunks generated" >> {log}
            exit 1
        fi
        echo "Filtered chunks generated successfully" >> {log}

        # Step 5: Concatenate the header and filtered chunks into a single SAM file
        cat {params.output_dir}/sam_header.sam {params.chunk_prefix}*.out > {params.filtered_sam}
        if [ $? -ne 0 ]; then
            echo "Error concatenating SAM files" >> {log}
            exit 1
        fi
        echo "Concatenated SAM files successfully" >> {log}

        # Step 6: Convert filtered SAM back to BAM
        samtools view -b {params.filtered_sam} > {output.filtered_bam} 2>> {log}
        if [ $? -ne 0 ]; then
            echo "Error in converting filtered SAM to BAM" >> {log}
            exit 1
        fi
        echo "Converted filtered SAM to BAM successfully" >> {log}

        # Clean up temporary files
        rm {params.chunk_prefix}*
        rm {params.filtered_sam} {params.inter_sam} {params.output_dir}/sam_header.sam
        echo "Cleaned up temporary files" >> {log}
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
        grep -v "JAFEMN" {params.output_dir}/{wildcards.pairs}.tmp.sam | grep -v "^@" >> {params.output_dir}/{wildcards.pairs}.header.sam
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





