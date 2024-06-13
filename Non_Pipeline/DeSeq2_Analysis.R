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
MetaData$Year <- as.character(MetaData$Year)


###############################
# Running DeSeq2 
###############################
#Makes DeSeq data set
dds <- DESeqDataSetFromMatrix(countData = Featurecount,
                              colData = MetaData,
                              design= ~ Treatment + Clone + Endophyte + Year + Month, tidy = TRUE)
#Run DeSeq function, 
dds_og <- DESeq(dds)
 dds <- dds_og

# filter for genes that have 10 occurrences in 1/4 the samples
keep <- rowSums(counts(dds) >= 1) >= (ncol(dds) / 4)
dds <- dds[keep, ]
res <- results(dds)

####################################
# Getting many comparisons
####################################

# Get Main Set of Results
res <- res[order(res$padj),]
summary(res)
head(res)

# Control vs Heat
Control_vs_Heat <- results(dds, contrast = c("Treatment", "Control", "Heat"))
Control_vs_Heat <- Control_vs_Heat[order(Control_vs_Heat$padj),]
summary(Control_vs_Heat)
head(Control_vs_Heat)

#
res <- res[order(res$padj),]
summary(res)
head(res)

#
res <- res[order(res$padj),]
summary(res)
head(res)

#
res <- res[order(res$padj),]
summary(res)
head(res)




# Scatter Plot
# par(mfrow=c(2,3))
# plotCounts(dds, gene="gene.117427.469.0", intgroup="Treatment")
# plotCounts(dds, gene="gene.94951.99.0", intgroup="Treatment")
# plotCounts(dds, gene="gene.122487.2597.0", intgroup="Treatment")
# plotCounts(dds, gene="gene.17870.375.2", intgroup="Treatment")
# plotCounts(dds, gene="gene.54571.320.0", intgroup="Treatment")
# plotCounts(dds, gene="gene.1127.808.0", intgroup="Treatment")

# Volcano Plot (template)
res$significance <- "Not Significant"
res$significance[res$padj < 0.05 & res$log2FoldChange > 2] <- "Significant Upregulated"
res$significance[res$padj < 0.05 & res$log2FoldChange < -2] <- "Significant Downregulated"

ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.8) + 
  scale_color_manual(values = c("Not Significant" = "grey", 
                                "Significant Upregulated" = "red", 
                                "Significant Downregulated" = "blue")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value") +
  theme_minimal() +
  theme(legend.position = "top")


# Volcano plot Control v Heat
Control_vs_Heat$significance <- "Not Significant"
Control_vs_Heat$significance[Control_vs_Heat$padj < 0.05 & Control_vs_Heat$log2FoldChange > 2] <- "Significant Upregulated"
Control_vs_Heat$significance[Control_vs_Heat$padj < 0.05 & Control_vs_Heat$log2FoldChange < -2] <- "Significant Downregulated"

ggplot(Control_vs_Heat, aes(x = log2FoldChange, y = -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.8) + 
  scale_color_manual(values = c("Not Significant" = "grey", 
                                "Significant Upregulated" = "red", 
                                "Significant Downregulated" = "blue")) +
  labs(title = "Control vs Heat", x = "Log2 Fold Change", y = "-Log10 P-value") +
  theme_minimal() +
  theme(legend.position = "top")












