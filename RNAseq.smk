# Author: Darrian Talamantes
# Affiliation: University of Georgia



# 1. Use fastqc and multiqc before running this. I could not get multiqc to work.
# 2. trim: This will use trimgalore to trim the reads.
# 3. make_transcriptome: Uses STAR to index and align reads to genome. It then makes lots of small bam files from the large one
# 4. scallop: Uses Scallop to create the transcriptome 
# 5. feature_counts: Uses featurecounts in R to do feature counting. 


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
align = config["directories"]["align"]
sep_bams = config["directories"]["sep_bams"]
# Files
genome = config["genome"]
star_manifest = config["star_manifest"]

# =================================================================================================
#   Functions
# =================================================================================================
def read_sample_names(file_path):
    with open(file_path, "r") as file:
        sample_names = [line.strip().strip('.') for line in file]
    return sample_names



# =================================================================================================
#   List Variables
# =================================================================================================

# Listing Progeny Files
sample_names_file = "/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/samples_list.txt"
SAMPLES = read_sample_names(sample_names_file)

paired_list_file = config["lists"]["paired_list_file"]
PAIRS = read_sample_names(paired_list_file)

# print("PAIRS:", PAIRS)  # Debugging: print the PAIRS list to verify contents

# =================================================================================================
#    Star Index files
# =================================================================================================

star_index_files = config["star_index_files"]

# =================================================================================================
#    Rules
# =================================================================================================


rule all:
    input:
        # expand(trimmed + "/{pairs}R1.fq.gz", pairs=PAIRS),
        # expand(trimmed + "/{pairs}R1_trimming_report.txt", pairs=PAIRS),
        # expand(trimmed + "/{pairs}R2.fq.gz", pairs=PAIRS),
        # expand(trimmed + "/{pairs}R2_trimming_report.txt", pairs=PAIRS)

        expand(config["kraken"]["classified"] + "/krakened_{pairs}R1.fq.gz", pairs=PAIRS)
        expand(config["kraken"]["classified"] + "/krakened_{pairs}R2.fq.gz", pairs=PAIRS)
        expand(config["kraken"]["unclassified"] + "/krakened_{pairs}R1.fq.gz", pairs=PAIRS)
        expand(config["kraken"]["unclassified"] + "/krakened_{pairs}R2.fq.gz", pairs=PAIRS) # Outputs for kraken

        # expand(config["directories"]["genome_idx"] + "/" + "{file}", file=star_index_files), # For indexing genome
        # config["directories"]["star_bams"] + "Aligned.sortedByCoord.out.bam", # mapping
        # config["directories"]["star_bams"] + "Log.out", # mapping
        # config["directories"]["star_bams"] + "Log.final.out", # mapping
        # config["directories"]["star_bams"] + "SJ.out.tab", # mapping
        # expand(config["directories"]["sep_bams"] + "{pairs}Aligned.sortedByCoord.out.bam", pairs=PAIRS), # mapping 2
        # gtf=config["scallop"]["output_file"] # Activates scallop2
        # counts = config["directories"]["features"] + "feature_counts.txt" # feature counts
        # pep_file = config["transdecoder"]["pep"], # predict makes this
        # fasta_gff3 = config["transdecoder"]["fasta_gff3"], # predict makes this     
        # pep_file_clean = config["transdecoder"]["pep_clean"]
   


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


## Look at the multiqc file and drop any that dont look good, then run the rest of the rules

include: "rules/trim.smk"
# # Here we run fastqc and multiqc manually. I will trim any samples with too many reads.

include: "rules/kraken.smk"
# include: "rules/star.smk"
# include: "rules/Fungal_removal.smk"
# include: "rules/scallop.smk"
# include: "rules/feature_counts.smk"
# include: "rules/annotation.smk"



