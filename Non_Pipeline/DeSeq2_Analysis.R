# This is R script used to analyze the genes within the RNA seq data
# Auther: Darrian Talamantes

# Notes:
# Before starting this script I deleted the suffix within the feature counts table

# Install Bioconductor if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install DESeq2 package from Bioconductor
BiocManager::install("DESeq2")



library(DESeq2)
library(ggplot2)

# File locations
MetaData_loc <- "/home/drt06/Documents/Tall_fescue/RNA_seq_fescue/Non_Pipeline/Meta_Data.csv"
Featurecount_loc <- "/home/drt06/Documents/Tall_fescue/Plant_Info/RNA_SEQ_backup_Data/feature_counts_namefix.txt"

###############################
# loading data
###############################
MetaData <- read.table(MetaData_loc, header = TRUE, sep = "," )
Featurecount <- read.table(Featurecount_loc, header = TRUE)

###############################
# Fixing data
###############################
Featurecount$feature <- rownames(Featurecount)
# Move the new column to the first position
Featurecount <- Featurecount[, c("feature", setdiff(names(Featurecount), "feature"))]


###############################
# Running DeSeq2 
###############################
#Makes DeSeq data set
dds <- DESeqDataSetFromMatrix(countData = Featurecount,
                              colData = MetaData,
                              design= ~ Treatment, tidy = TRUE )
#Run DeSeq function
dds <- DESeq(dds)

#Look at results table
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res)



