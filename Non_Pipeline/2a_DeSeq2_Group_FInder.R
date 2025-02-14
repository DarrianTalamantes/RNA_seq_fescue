#Purpose: This script will create files that has the upregulated and 
# Downregulated genes from our different treatments.
# These files will then be put into GoTerm analysis




#libraries
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
library(unix)

# Increase memory of R to 12 GB
rlimit_as(1e12)



###############################
# loading data
###############################
PlantInfo="/home/darrian/Documents/Plant_Info"
RNA_seq_fescue="/home/darrian/Documents/UGA/RNA_seq_fescue"

MetaData_loc <- paste0(RNA_seq_fescue, "/Non_Pipeline/Meta_Data.csv")
Featurecount_loc <- paste0(PlantInfo, "/RNA_SEQ_backup_Data/feature_counts_namefix.txt")


MetaData <- read.table(MetaData_loc, header = TRUE, sep = "," )
Featurecount <- read.table(Featurecount_loc, header = TRUE)



