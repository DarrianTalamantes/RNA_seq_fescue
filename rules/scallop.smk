# This rule is to mainly use scallop2 however before I use it I 
# recreate the header, index the bam file and split it by chromosome.
# I then can speed up scallop by running it on smaller bam files.


# Recreates the header and adds it to the .bam file
rule fix_bam_header:
    input:
        bam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.bam",
        ref = config["genome"]
    output:
        bam_fixed = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered_fixed.out.bam"
    conda:
        "../Conda_Envs/smatools.yaml"
    log:
        "logs/fix_bam_header.log"
    shell:
        """
        samtools faidx {input.ref}
        cut -f1,2 {input.ref}.fai | awk '{{print "@SQ\\tSN:"$1"\\tLN:"$2}}' > new_header.sam
        samtools reheader new_header.sam {input.bam} > {output.bam_fixed} 2> {log}
        """

# Indexes the new bam file
rule index_bam:
    input:
        bam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered_fixed.out.bam"
    output:
        bai = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered_fixed.out.bam.bai"
    conda:
        "../Conda_Envs/smatools.yaml"
    log:
        "logs/index_bam.log"
    shell:
        """
        samtools index {input.bam} 2> {log}
        """

# Splits the bam file by chromosome
checkpoint split_bam_by_chr:
    input:
        bam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered_fixed.out.bam",
        bai = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered_fixed.out.bam.bai"
    output:
        chrom = config["directories"]["big_bam_chrom"]
    conda:
        "../Conda_Envs/smatools.yaml"
    log:
        "logs/split_bam_by_chr.log"
    shell:
        """
        mkdir -p {output}
        samtools idxstats {input.bam} | cut -f1 | grep -v '*' | while read chr; do
            samtools view -b {input.bam} $chr > {output.chrom}/${{chr}}.bam
        done 2> {log}
        """
def get_bam_files(wildcards):
    checkpoint_output = checkpoints.split_bam_by_chr.get(**wildcards).output[0]
    bam_files = glob.glob(f"{checkpoint_output}/*.bam")
    return bam_files



# Runs scallop2 on all chromosome files
rule scallop2:
    input:
        bam = get_bam_files
    conda:
        "../Conda_Envs/scallop2.yaml"
    threads: config["scallop"]["threads"]
    output:
        gtf = config["directories"]["scallop_out"] + "/{chrom}.gtf"
    log:
        "logs/scallop2_{chrom}.log"
    shell:
        """
        mkdir -p scallop2_output
        scallop2 --num-threads {threads} -i {input.bam} -o {output.gtf} 2> {log}
        """

# merges the gtf chromosome files into one file, this may produce a huuuge log file. Its not helpfull other than knowing scallop is running
rule merge_gtfs:
    input:
        gtfs = lambda wildcards: glob.glob(config["directories"]["scallop_out"] + "/*.gtf")
    output:
         gtf = config["scallop"]["output_file"]
    conda:
       "../Conda_Envs/gtftools.yaml"
    shell:
        """
        cat {input.gtfs} | grep -v '^#' | sort -k1,1 -k4,4n > {output.gtf}
        """


########################################################
# Old rule: This works but takes suuuper long.
########################################################
# rule scallop2:
#     input:
#         bam = config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.bam"
#     conda:
#         "../Conda_Envs/scallop2.yaml"
#     log:
#         "logs/scallop2_{chrom}.log"
#     params:
#         threads = config["scallop"]["threads"]
#     threads: config["scallop"]["threads"]
#     output:
#         gtf = config["scallop"]["output_file"]
#     shell:
#         """
#         scallop2 --num-threads {params.threads} -i {input.bam} -o {output.gtf} 2> {log}
#         """
