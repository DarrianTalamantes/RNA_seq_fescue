# Purpose: This is to be ran after bam_read_counter.sh. This counts the reads that aligned
# to the epichloe genome and this R script graphs it by E+ or E-.
# Auther: Darrian Talamantes

# Notes:
# Before starting this script I deleted the suffix within the feature counts table

# Install Bioconductor if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(grid)
library(data.table)
library(pheatmap)
library(variancePartition)
library(BiocGenerics) 
library(lme4)        
library(gridExtra)
library(patchwork)
library(cowplot)



# install.packages("unix")

# Increase memory of R to 12 GB
memory.limit(size=12000)  #Might not need. Memory should not be capped in newer R versions

# File locations
data_folder <- "/home/darrian/Documents/RNA_seq_fescue/Lists"
bam_counts_epi_loc <- paste0(data_folder, "/bam_read_counts_clean.txt")
bam_counts_fesc_loc <- paste0(data_folder, "/bam_read_counts_fescue_clean.txt")
bam_counts_all_loc <- paste0(data_folder, "/bam_read_counts_fesc_epi_clean.txt")

###############################
# loading data
###############################
bam_counts_epi <- read.table(bam_counts_epi_loc, header = FALSE)
bam_counts_fesc <- read.table(bam_counts_fesc_loc, header = FALSE)
bam_counts_all <- read.table(bam_counts_all_loc, header = FALSE)



################################################################################
# Fixing data and making graph of EPichloe neg to pos
################################################################################


bam_counts_epi <- bam_counts_epi %>%
  mutate(
    # Capture CTE## and N/P together
    matches = str_match(V1, "(CTE[0-9]+)([NP])"),
    Genotype = matches[,2],  # CTE##
    Epichloe = matches[,3]    # N or P
  ) %>%
  select(-matches)




ggplot(bam_counts_epi, aes(x = Genotype, y = V2, color = Epichloe)) +
  geom_jitter(width = 0.2, height = 0, size = 3) +  # jitter points slightly on X-axis
  theme_bw() +
  labs(x = "Genotype", y = "Read Count", color = "Epichloe") +
  ggtitle("Epichloe Read Counts Per Sample") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18), 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  scale_color_manual(values = c("N" = "#1b9e77", "P" = "#d95f02"))

################################################################################
# Comparing epichloe reads to fescue reads and all reads
################################################################################

bam_counts_epi_colnames <- c("Sample", "Epichloe Reads")
bam_counts_fesc_colnames <- c("Sample", "Fescue Reads")
bam_counts_all_colnames <- c("Sample", "Total Reads")

colnames(bam_counts_epi) <- bam_counts_epi_colnames
colnames(bam_counts_fesc) <- bam_counts_fesc_colnames
colnames(bam_counts_all) <- bam_counts_all_colnames

df_list <- list(bam_counts_epi, bam_counts_fesc, bam_counts_all)
merged <- Reduce(function(x, y) merge(x, y, by = "Sample", all = TRUE), df_list)
merged$`Percent Epichloe` <- round(merged$`Epichloe Reads`/merged$`Total Reads`,5)
merged$`Percent Fescue` <- round(merged$`Fescue Reads`/merged$`Total Reads`,5)

write.csv(merged, paste0(data_folder, "/Read_count_by_organism.csv"), row.names = FALSE)


