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
(SAMPLES,) = glob_wildcards(datadir + "{sample}.fastq.gz")





# =================================================================================================
#    Rules
# =================================================================================================


rule all:
    input:
        "multiqc_report.html"

# This rule runs fastqc on all data fastq files
rule fastqc:
    input:
        fastq = expand("{datadir}/{sample}.fastq.gz", sample = SAMPLES)
    output:
        html = qcdir + "/{sample}_fastqc.html",
        zip = qcdir + "/{sample}_fastqc.zip"
    log:
         f"{logs}/{sample}_fastqc.log"
    conda:
        "Conda_Env/multiqc.yaml"
    shell:
        """
        fastqc  -o {qcdir} {input.fastq} 
        echo "fastqc  -o {qcdir} {input.fastq}" >> {log}
        """

# This rule uses multiqc on the output of the rul fastqc
rule multiqc:
    input:
        expand(qcdir + "/{sample}_fastqc.html")
    output:
        "multiqc_report.html"
    conda:
        "Conda_Env/multiqc.yaml"        
    shell:
        "multiqc {qcdir} -o {qcdir} --filename multiqc_report.html"




        