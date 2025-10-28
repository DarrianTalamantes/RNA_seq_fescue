# Author: Darrian Talamantes
# Affiliation: University of Georgia

# Purpose of program: This will use the *Aligned.sortedByCoord.out.bam files to filter out the epichloe reads. 
# From there I will make a VCF file our of them and analyse them in R

# =========================================================================================================
#   Importing wrapper stuff
# =========================================================================================================
from snakemake.io import expand
import glob


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

for path in config["directories"].values():
    os.makedirs(path, exist_ok=True)


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
star_index_files = list(star_index_files)  # Ensure it's a list


# =================================================================================================
#    Rules
# =================================================================================================
print(expand(config["directories"]["genome_idx"] + "/{file}", file=star_index_files))

print(config["directories"]["big_bam"])
print(config["directories"]["filtered_bam_big"])
print(config["fungal_removal"]["threads"])

# You can not have anymore that one commented out line when defining your inputs here, There can not be blank line after input.
rule all:
    input:
        # # Fungal Filtering
        config["directories"]["filtered_bam_big"] + "/Aligned.sortedByCoord_filtered.out.bam",
        expand(config["directories"]["filtered_bams"] + "/{pairs}Aligned.sortedByCoord_filtered.out.bam", pairs=PAIRS)

        # # # Scallop and feature counts
        # config["scallop"]["output_file_big"]

        # # # Feature counts
        # config["directories"]["features"] + "feature_counts.txt"
        



   



## Look at the multiqc file and drop any that dont look good, then run the rest of the rules

# This is the order in which we use the programs.


include: "rules/fungal_removal.smk" # Dry runs may show that concatenate_and_convert_big will not work. It will

# include: "rules/scallop.smk"

# include: "rules/feature_counts.smk"



