# Author: Darrian Talamantes
# Affiliation: University of Georgia


# To use this id suggest doing one step at a time. The rule all has code seperated into blocks and below it you can include the rules
# you want to use

# 1. Use fastqc and multiqc before running this. I could not get multiqc to work.
# 2. trim: This will use trimgalore to trim the reads.
# 3. make_transcriptome: Uses STAR to index and align reads to genome. It then makes lots of small bam files from the large one
# 4. scallop: Uses Scallop to create the transcriptome 
# 5. feature_counts: Uses featurecounts in R to do feature counting. 


# =========================================================================================================
#   Importing wrapper stuff
# =========================================================================================================
from snakemake.io import expand


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

#This is the base name of all the files without R1 and R2
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

# You can not have anymore that one commented out line when defining your inputs here
rule all:
    input:

        # # Kraken outputs
        expand(config["kraken"]["classified"] + "/krakened_{pairs}.txt", pairs=PAIRS),
        # expand(config["kraken"]["fungal"] + "/{pairs}R1.fq", pairs=PAIRS),
        # expand(config["kraken"]["fungal"] + "/{pairs}R2.fq", pairs=PAIRS),
        # expand(config["kraken"]["non_fungal"] + "/{pairs}R1.fq", pairs=PAIRS),
        # expand(config["kraken"]["non_fungal"] + "/{pairs}R2.fq", pairs=PAIRS)
        # config["kraken"]["db_name"] + "/hash.k2d",   # This line is probs useless if you download the db     


        # # Star big bam
        # expand(config["directories"]["genome_idx"] + "/{file}", file=star_index_files),
        # config["directories"]["big_bam"] + "Aligned.sortedByCoord.out.bam",
        # config["directories"]["big_bam"] + "Log.out",
        # config["directories"]["big_bam"] + "Log.final.out",
        # config["directories"]["big_bam"] + "SJ.out.tab",

        # # Fungal Filtering
        # expand(config["directories"]["filtered_bams"] + "/{pairs}Aligned.sortedByCoord_filtered.out.bam", pairs=PAIRS)
        # config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.bam"

        

        # # Star sep bams
        # expand(config["directories"]["sep_bams"] + "{pairs}Aligned.sortedByCoord.out.bam", pairs=PAIRS),
        # expand(config["directories"]["sep_bams"] + "{pairs}Log.out", pairs=PAIRS),
        # expand(config["directories"]["sep_bams"] + "{pairs}Log.final.out", pairs=PAIRS),
        # expand(config["directories"]["sep_bams"] + "{pairs}SJ.out.tab", pairs=PAIRS),

        # # Scallop
        # config["scallop"]["output_file"],
        # config["directories"]["features"] + "feature_counts.txt"


        # Stuff I have not totally labeled yet
        # gtf=config["scallop"]["output_file"] # Activates scallop2
        # counts = config["directories"]["features"] + "feature_counts.txt" # feature counts
        # pep_file = config["transdecoder"]["pep"], # predict makes this
        # fasta_gff3 = config["transdecoder"]["fasta_gff3"], # predict makes this     
        # pep_file_clean = config["transdecoder"]["pep_clean"]

        # # Trimming   
        # expand(trimmed + "/{pairs}R1.fq.gz", pairs=PAIRS), 
        # expand(trimmed + "/{pairs}R1_trimming_report.txt", pairs=PAIRS),
        # expand(trimmed + "/{pairs}R2.fq.gz", pairs=PAIRS), 
        # expand(trimmed + "/{pairs}R2_trimming_report.txt", pairs=PAIRS)


   



## Look at the multiqc file and drop any that dont look good, then run the rest of the rules

# This is the order in which we use the programs.
# include: "rules/trim.smk"
# # Here we run fastqc and multiqc manually. I will trim any samples with too many reads by just cutting them to a length of the next largest file
# # The file name is RunFastQC.sh then multiqc to see what needs to be cut down to length 

include: "rules/kraken.smk"
# include: "rules/star.smk"
# include: "rules/fungal_removal.smk"

# include: "rules/scallop.smk"
# include: "rules/feature_counts.smk"

# include: "rules/annotation.smk"



