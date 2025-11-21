Steps to rerun Scripts.

#############################################################################

# Notes for running analysis when Darrian leaves

1. First thing you must do is run the snakemake pipeline. I never got it to work on its own. In the RNAseq.smk there is a list of rules that you must run in order. Other things you will need to make sure the pipeline runs is a combined epichloe and tall fescue genome and all the reads from the rna seq data. Necessary locations for these can be found in the config.yaml

Below is the rules in order they need to be ran as a backup

include: "rules/kraken.smk" # this rule only works with snakemake version: snakemake/6.9.1-Mamba-4.11.0-4
include: "rules/star.smk" # Switch to snakemake version: snakemake/7.22.0-foss-2022a
include: "rules/fungal_removal.smk" # Dry runs may show that concatenate_and_convert_big will not work. It will
include: "rules/scallop.smk"
include: "rules/feature_counts.smk"
include: "rules/annotation.smk"

**Update**
I have made this into individual sapelo2 submissions. There is now a folder called 1_RNA_SMK_submissions. This should allow you to simply submit the RNA seek pipline as individual scripts to sapelo2. You should be able to submit 2 and 3 at the same time but all others must be completed in order. 

I forget the exact amount of memory everything used but I think what I set stuff to should be fine. If I remmeber correctly it was creating the big bam file that used the most amount of memory.


2. Once you run the snakemake pipeline you will have two main output files.
- Fescue_transcriptome.gtf
- feature_counts.txt

These output files are all you will need to run R files in the folder 2_Non_Pipeline. In order to move on to Go-term analysis with the EnTAP pipeline you must run up to 2d_Upset_plots.R

3. For the EnTAP pipline ran in 3a_GoaTools.py you need the output of 2d_Upset_plots.R. The ouputs of this look like
- Treatments_Up_Down_reg.csv
- Genotypes_Up_Down_reg_HeatvsControl.csv 
- ect.

I never automated this 3a_GoaTools.py script. Thus you must run it and change the line

DEG_count_Data = "/home/darrian/Documents/RNA_seq_fescue/r_data/Treatments_Up_Down_reg.csv" # Switch this line to change DEG set

in order to analyze the set of DEGs that you care about.

4. The directory 4_Graph_making makes a few of the graphs for publication. Mainly better heatmaps and upset plots.

5. The directory 5_Epichloe_checker is full of scripts that will check the read counts of bam files that align to tall fescue, epichloe, and the two genomes combined. I then have an R script in there that. To run this you have to run the EpiSnake.smk pipeline first to remove epichloe reads as the origional pieline does not save these reads form the combined bam files. After this run 5a-5c which will count the reads.