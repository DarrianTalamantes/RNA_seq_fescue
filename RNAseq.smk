# Author: Darrian Talamantes
# Affiliation: University of Georgia



# 1. QC/QA using mutliQC and Trimmomatic to remove adapters and low quality reads 
# 2. Map to genome using Salmon (quasi mapper) it also quantifies the reads 
# 3. Differential analysis using DESeq2 (via R) 

# =========================================================================================================
#     Load config file
# =========================================================================================================
configfile: "config.yaml"

# Config importing paths to directories 
datadir = config["directories"]["datadir"]
qcdir = config["directories"]["qcdir"]
logs = config["directories"]["logs"]

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




# =================================================================================================
#    Rules
# =================================================================================================


rule all:
    input:
        qcdir + "multiqc_report.html"

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

# This rule uses multiqc on the output of the rule fastqc
# rule multiqc:
#     input:
#         expand("{qcdir}/{sample}_fastqc.html", qcdir=qcdir, sample=SAMPLES)
#     output:
#         "multiqc_report.html"
#     conda:
#         "Conda_Envs/multiqc.yaml"        
#     shell:
#         "multiqc {qcdir} -o {qcdir} --filename multiqc_report.html"


rule multiqc:
    input:
        expand(qcdir + "/{sample}_fastqc.html", sample=SAMPLES)
    output:
        html_report = qcdir + "/multiqc_report.html"
    params:
        outdir = qcdir,
        options = "--filename multiqc_report.html"
    message: 
        "Performing MultiQC on the FastQC results"
    wrapper:
        "v1.20.0/bio/multiqc" 
               