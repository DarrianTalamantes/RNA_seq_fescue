# Objective: This rule should go between STAR and Scallop. Once we have mapped everything to both genomes this will seperate fungal reads.

##################################################
# Getting rid of fungus
##################################################

rule split_and_filter_big:
    input:
        big_bam = config["directories"]["big_bam"] + "Aligned.sortedByCoord.out.bam"
    conda:
        "../Conda_Envs/samtools.yaml"
    params:
        output_dir = config["directories"]["filtered_bam_big"],
        chunk_prefix = config["directories"]["filtered_bam_big"] + "/chunk_",
        inter_sam = config["directories"]["big_bam"] + "Aligned.sortedByCoord.out.sam",
        lines_per_chunk = config["fungal_removal"]["lines_per_chunk"]
    threads: 32
    output:
        chunked_outs = directory(config["directories"]["filtered_bam_big"]),  # Use directory output
        header = config["directories"]["filtered_bam_big"] + "/sam_header.sam"  # Ensure header is created
    log:
        "logs/split_and_filter_big.log"
    shell:
        """
        echo "Starting split_and_filter_big rule" >> {log}

        # Step 1: Convert BAM to SAM
        samtools view -h {input.big_bam} > {params.inter_sam} 2>> {log}
        echo "Converted BAM to SAM successfully" >> {log}

        # Step 2: Extract the SAM header
        grep "^@" {params.inter_sam} > {params.output_dir}/sam_header.sam 2>> {log}
        echo "Extracted SAM header successfully" >> {log}

        # Step 3: Split the SAM file into chunks
        grep -v "^@" {params.inter_sam} > temp_sam_no_header.sam
        split -l {params.lines_per_chunk} temp_sam_no_header.sam {params.chunk_prefix} 2>> {log}
        echo "Split SAM file into chunks successfully" >> {log}

        # Step 4: Use parallel to grep and filter each chunk
        ls {params.chunk_prefix}* | parallel -j {threads} "grep -v 'JAFEMN' {{}} > {{}}.out" 2>> {log}
        echo "Filtered chunks successfully" >> {log}

        # Clean up intermediate files
        rm temp_sam_no_header.sam {params.inter_sam}
        echo "Cleaned up intermediate files" >> {log}
        """

wildcards_dict = glob_wildcards(config["directories"]["filtered_bam_big"] + "/chunk_{i}.out")
chunk_count = len(wildcards_dict.i)

rule concatenate_and_convert_big:
    input:
        filtered_chunks=lambda wildcards: expand(
            config["directories"]["filtered_bam_big"] + "/chunk_{i}.out", 
            i=glob_wildcards(config["directories"]["filtered_bam_big"] + "/chunk_{i}.out").i
        ),
        header=lambda wildcards: config["directories"]["filtered_bam_big"] + "/sam_header.sam"
    conda:
        "../Conda_Envs/samtools.yaml"
    params:
        filtered_sam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.sam"
    output:
        filtered_bam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.bam"
    log:
        "logs/concatenate_and_convert_big.log"
    shell:
        """
        echo "Starting concatenate_and_convert_big rule" >> {log}

        # Step 5: Concatenate filtered chunks
        echo "Concatenating filtered chunks..." >> {log}
        cat {input.header} {input.filtered_chunks} > {params.filtered_sam} 2>> {log}
        if [ $? -ne 0 ]; then
            echo "Error during concatenation" >> {log}
            exit 1
        fi
        echo "Concatenated filtered chunks successfully" >> {log}

        # Step 6: Convert filtered SAM to BAM
        samtools view -b {params.filtered_sam} > {output.filtered_bam} 2>> {log}
        echo "Converted filtered SAM to BAM successfully" >> {log}

        # Step 7: Clean up intermediate files
        rm {params.filtered_sam}
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





