# Author: Darrian Talamantes
# Affiliation: University of Georgia



# 1. QC/QA using mutliQC and Trim Galore to remove adapters and low quality reads 
# 2. Create a transcriptome using trinity from both genome and sequence data
# 3. Use salmon to quantify shit 
# 4. Differential analysis using DESeq2 (via R) 


# =========================================================================================================
#   Importing wrapper stuff
# =========================================================================================================


# =========================================================================================================
#     Load config file
# =========================================================================================================
configfile: "config.yaml"

# Config importing paths to directories 
datadir = config["directories"]["datadir"]
qcdir = config["directories"]["qcdir"]
logs = config["directories"]["logs"]
mqcdir = config["directories"]["mqcdir"]
trimmed = config["directories"]["trimmed"]
genome_idx = config["directories"]["genome_idx"]
align = config["directories"]["Alignment"]

# Files
genome = config["genome"]


# =================================================================================================
#   Functions
# =================================================================================================
def read_sample_names(file_path):
    with open(file_path, "r") as file:
        sample_names = [line.strip() for line in file]
    return sample_names



# =================================================================================================
#   List Variables
# =================================================================================================

# Listing Progeny Files
sample_names_file = "/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/samples_list.txt"
SAMPLES = read_sample_names(sample_names_file)

paired_list_file = "/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/paired_list.txt"
PAIRS = read_sample_names(paired_list_file)


# =================================================================================================
#    Rules
# =================================================================================================


rule all:
    input:
        # mqcdir + "/multiqc_report.html"
        fasta_fwd= expand(trimmed + "/{pairs}R1.fq.gz", pairs=PAIRS),
        report_fwd= expand(trimmed + "/{pairs}R1_trimming_report.txt", pairs=PAIRS),
        fasta_rev= expand(trimmed + "/{pairs}R2.fq.gz", pairs=PAIRS),
        report_rev= expand(trimmed + "/{pairs}R2_trimming_report.txt", pairs=PAIRS),
        directory(config["directories"]["genome_idx"])


# This rule runs fastqc on all data fastq files
rule fastqc:
    input:
        fastq = datadir + "/{sample}.fastq.gz"
    output:
        html = qcdir + "/{sample}_fastqc.html",
        zip = qcdir + "/{sample}_fastqc.zip"   
    conda:
        "Conda_Envs/multiqc.yaml"
    shell:
        """
        fastqc  -o {qcdir} {input.fastq} 
        """

#################################
# Multiqc (could not get multiqc installed in conda)
#################################
# rule multiqc:
#     input:
#         expand(qcdir + "/{sample}_fastqc.html", sample=SAMPLES),
#     output:
#         mqcdir + "/multiqc_report.html",
#         mqcdir + "/multiqc_report.zip",
#     params:
#         extra= "--verbose",
#     conda:
#         "Conda_Envs/multiqc.yaml"
#     wrapper:
#         "v3.9.0/bio/multiqc"
####
# Code below was ran instead of the multiqc rule
# multiqc . 
###

include: "rules/trim.smk"
include: "make_transcriptome.smk"




