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
#   List Variables
# =================================================================================================

# Listing Progeny Files
sample_names_file = "sample_list.txt"
SAMPLES = read_sample_names(sample_names_file)




# =================================================================================================
#    Rules
# =================================================================================================


rule all:
    input:
        "multiqc_report.html"

# This rule runs fastqc on all data fastq files
rule fastqc:
    input:
        fastq = datadir + "/{sample}.fastq.gz"
    output:
        html = qcdir + "/{sample}_fastqc.html",
        zip = qcdir + "/{sample}_fastqc.zip"   
    conda:
        "Conda_Env/multiqc.yaml"
    shell:
        """
        fastqc  -o {qcdir} {input.fastq} 
        """

# This rule uses multiqc on the output of the rule fastqc
rule multiqc:
    input:
        expand("{qcdir}/{sample}_fastqc.html", qcdir=qcdir, sample=SAMPLES)
    output:
        "multiqc_report.html"
    conda:
        "Conda_Env/multiqc.yaml"        
    shell:
        "multiqc {qcdir} -o {qcdir} --filename multiqc_report.html"


def read_sample_names(file_path):
    with open(file_path, "r") as file:
        sample_names = [line.strip() for line in file]
    return sample_names

        