# Author: Darrian Talamantes
# This script will load in goaterm files and make charts of them


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
library(ComplexUpset)

# File locations
data_folder <- "/home/darrian/Documents/RNA_seq_fescue/Goatools_data"

# E+ and E- comparison when looking at genotypes
gene_counts = read.csv(paste0(data_folder, "/CTE31_Down_results.txt"), header = TRUE) # from 3a

# E+ and E- comparison when looking at treatment groups
gene_counts = read.csv(paste0(data_folder, "/Control_Down_results.txt"), header = TRUE) # from 3a



