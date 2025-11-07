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
bam_counts_loc <- paste0(data_folder, "/bam_read_counts.txt")

###############################
# loading data
###############################
bam_counts <- read.table(bam_counts_loc, header = FALSE)



################################################################################
# Fixing data and making graph
################################################################################


bam_counts <- bam_counts %>%
  mutate(
    # Capture CTE## and N/P together
    matches = str_match(V1, "(CTE[0-9]+)([NP])"),
    Genotype = matches[,2],  # CTE##
    Epichloe = matches[,3]    # N or P
  ) %>%
  select(-matches)




ggplot(bam_counts, aes(x = Genotype, y = V2, color = Epichloe)) +
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





