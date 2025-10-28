# This rule is to mainly use scallop2 however before I use it I 
# recreate the header, index the bam file and split it by chromosome.
# I then can speed up scallop by running it on smaller bam files.


# # Recreates the header and adds it to the .bam file
# rule fix_bam_header:
#     input:
#         bam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.bam",
#         ref = config["genome"]
#     output:
#         bam_fixed = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered_fixed.out.bam"
#     conda:
#         "../Conda_Envs/samtools.yaml"
#     log:
#         "logs/fix_bam_header.log"
#     shell:
#         """
#         samtools faidx {input.ref}
#         cut -f1,2 {input.ref}.fai | awk '{{print "@SQ\\tSN:"$1"\\tLN:"$2}}' > new_header.sam
#         samtools reheader new_header.sam {input.bam} > {output.bam_fixed} 2> {log}
#         """

# # Indexes the new bam file
# rule index_bam:
#     input:
#         bam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered_fixed.out.bam"
#     output:
#         bai = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered_fixed.out.bam.bai"
#     conda:
#         "../Conda_Envs/samtools.yaml"
#     log:
#         "logs/index_bam.log"
#     shell:
#         """
#         samtools index -c {input.bam} 2> {log}
#         """

# # Splits the bam file by chromosome
# checkpoint split_bam_by_chr:
#     input:
#         bam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered_fixed.out.bam",
#         bai = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered_fixed.out.bam.bai"
#     output:
#         chrom = directory(config["directories"]["big_bam_chrom"])
#     conda:
#         "../Conda_Envs/samtools.yaml"
#     log:
#         "logs/split_bam_by_chr.log"
#     shell:
#         """
#         mkdir -p {output.chrom}
#         samtools idxstats {input.bam} | cut -f1 | grep -v '*' | while read chr; do
#             # Split the BAM file by chromosome, keep header, and sort by reference position (-o)
#             samtools view -b -h {input.bam} $chr | samtools sort -o {output.chrom}/${{chr}}.bam
#             samtools index -c {output.chrom}/${{chr}}.bam
#         done 2> {log}
#         """

# def get_chromosomes(wildcards):
#     checkpoint_data = checkpoints.split_bam_by_chr.get(**wildcards)  # Wait for checkpoint
#     checkpoint_output = checkpoint_data.output.chrom
#     bam_files = glob.glob(f"{checkpoint_output}/*.bam")
#     return [os.path.basename(f).replace(".bam", "") for f in bam_files]


# # Get chromosome names from the split BAM checkpoint
# def get_chromosomes(wildcards):
#     checkpoint_data = checkpoints.split_bam_by_chr.get(**wildcards)  # Ensure BAMs are split
#     checkpoint_output = checkpoint_data.output.chrom  # Get directory with split BAMs
#     bam_files = glob.glob(f"{checkpoint_output}/*.bam")  # List BAM files
#     return [os.path.basename(f).replace(".bam", "") for f in bam_files]  # Extract chromosome names

# # Rule for transcript assembly with Scallop2
# rule scallop2:
#     input:
#         bam = lambda wildcards: f"{config['directories']['big_bam_chrom']}/{wildcards.chrom}.bam"
#     output:
#         gtf = config["directories"]["scallop_out"] + "/{chrom}.gtf"  
#     conda:
#         "../Conda_Envs/scallop2.yaml"
#     params:
#         directory = config["directories"]["scallop_out"]  
#     threads: config["scallop"]["threads"]
#     log:
#         "logs/scallop2_{chrom}.log"
#     shell:
#         """
#         mkdir -p {params.directory}
#         scallop2 --num-threads {threads} -i {input.bam} -o {output.gtf} 2> {log}
#         """

# # **Make merge_gtfs explicitly depend on scallop2 output**
# rule merge_gtfs:
#     input:
#         gtfs = expand(config["directories"]["scallop_out"] + "/{chrom}.gtf", chrom=get_chromosomes)
#     output:
#         gtf = config["scallop"]["output_file"]
#     shell:
#         """
#         cat {input.gtfs} | grep -v '^#' | sort -k1,1 -k4,4n > {output.gtf}
#         """

########################################################
# Old rule: This works but takes suuuper long.
########################################################
rule sort_bam_by_coord:
    input:
        bam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.bam"
    output:
        bam_sorted = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered_sorted.out.bam"
    conda:
        "../../Conda_Envs/samtools.yaml"
    params:
        backups = config["directories"]["bam_backup"] + "/Aligned.sortedByCoord_filtered_sorted.out.bam"
    log:
        "logs/sort_bam_by_coord_big.log"
    threads: config["scallop"]["threads"]
    shell:
        """
        samtools sort -@ {threads} -o {output.bam_sorted} {input.bam} 2> {log}
        cp {output.bam_sorted} {params.backups}
        """

rule index_bam:
    input:
        bam_sorted = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered_sorted.out.bam"
    output:
        csi = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered_sorted.out.bam.csi"
    conda:
        "../../Conda_Envs/samtools.yaml"
    log:
        "logs/index_bam_big.log"
    threads: config["scallop"]["threads"]
    shell:
        """
        samtools index -c {input.bam_sorted} 2> {log}
        """


rule scallop2_big:
    input:
        bam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered_sorted.out.bam",
        csi = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered_sorted.out.bam.csi"
    output:
        gtf = config["scallop"]["output_file_big"]
    conda:
        "../Conda_Envs/scallop2.yaml"
    log:
        "logs/scallop2_Big.log"
    threads: config["scallop"]["threads"]
    shell:
        """
        scallop2 --num-threads {threads} -i {input.bam} -o {output.gtf} > {log}
        """






# Change this Aligned.sortedByCoord_filtered to test_out for testing.