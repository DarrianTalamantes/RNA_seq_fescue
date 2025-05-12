# This is R script used to analyze the genes within the RNA seq data
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
# install.packages("unix")

# Increase memory of R to 12 GB
memory.limit(size=12000)  #Might not need. Memory should not be capped in newer R versions

# File locations
data_folder <- "/home/darrian/Documents/RNA_seq_fescue/r_data"
Featurecount_loc <- paste0(data_folder, "/feature_counts_name_fixed.txt")

###############################
# loading data
###############################
Featurecount <- read.table(Featurecount_loc, header = TRUE)

###############################
# Fixing data
###############################
# Featurecount$feature <- rownames(Featurecount)
# # Move the new column to the first position
# Featurecount <- Featurecount[, c("feature", setdiff(names(Featurecount), "feature"))]

# Creating the MetaData
IDs <- colnames(Featurecount)
IDs <- as.data.frame(IDs)
IDs <- IDs[-1, ]
# write.csv(IDs,"featurenames.txt",row.names = FALSE)


metadata <- data.frame(SampleName = IDs) %>%
  mutate(
    Clone = str_extract(SampleName, "CTE\\d+"),
    Endophyte = str_extract(SampleName, "CTE\\d+(N|P)"),
    Endophyte = str_sub(Endophyte, -1),
    
    Treatment = str_extract(SampleName, "_(H|HP|C)\\d"),
    Treatment = str_remove_all(Treatment, "_|\\d"),
    
    Replicate = str_extract(SampleName, "(H|HP|C)(\\d)"),
    Replicate = str_extract(Replicate, "\\d"),
    
    HarvestCode = str_extract(SampleName, "(O17|J16|O16|J17)"),
    
    Month = case_when(
      HarvestCode %in% c("J16", "J17") ~ "June",
      HarvestCode %in% c("O16", "O17") ~ "October",
      TRUE ~ NA_character_
    ),
    
    Year = case_when(
      HarvestCode == "J16" ~ "2016",
      HarvestCode == "O16" ~ "2016",
      HarvestCode == "J17" ~ "2017",
      HarvestCode == "O17" ~ "2017",
      TRUE ~ NA_character_
    ),
    
    HarvestTime = paste(Month, Year, sep = "_"),
    
    SampleID = str_extract(SampleName, "A\\d+")
  )

rep_str = c('N'='Negative','P'='Positive')
metadata$Endophyte <- str_replace_all(metadata$Endophyte, rep_str)

rep_str2 = c('HP'='HeatxPercipitation','H'='Heat', 'C' = 'Control')
metadata$Treatment <- str_replace_all(metadata$Treatment, rep_str2)
rep_str3 = c('Heateat'='Heat')
metadata$Treatment <- str_replace_all(metadata$Treatment, rep_str3)
metadata <- subset(metadata, select = -c(HarvestCode))


Featurecount <- Featurecount[, colnames(Featurecount) %in% metadata$SampleName]

ncol(Featurecount)
nrow(metadata)
###############################
# Running DeSeq2 
###############################
#Makes DeSeq data set
dds <- DESeqDataSetFromMatrix(countData = Featurecount,
                              colData = metadata,
                              design= ~ Clone + Month + Year + Endophyte + Treatment)

#Run DeSeq function, 
dds_og <- DESeq(dds)
 dds <- dds_og

 
# filter for genes that have 10 occurrences in 1/4 the samples
keep <- rowSums(counts(dds) >= 5) >= (ncol(dds) / 4)
dds <- dds[keep, ]

####################################
# Function to get results.
####################################
get_results <- function(DeSeqData,column,t1,t2){
  if(!is.character(column) || !is.character(t1) || !is.character(t2)) {
    stop("The column and treatment levels (t1, t2) must be provided as strings.")
  }
  output <- results(DeSeqData, contrast = c(column, t1, t2))
  output <- output[order(output$padj),]
  return(output)
}

####################################
# Getting many comparisons
####################################
# Control vs Heat
Control_vs_Heat <- get_results(dds,"Treatment","Control","Heat")
summary(Control_vs_Heat)
Control_vs_Heat <- Control_vs_Heat[order(Control_vs_Heat$pvalue),]
head(Control_vs_Heat)

# Control vs HeatxPercipitation
Control_vs_HxP <- get_results(dds,"Treatment","HeatxPercipitation","Control")
summary(Control_vs_HxP)
Control_vs_HxP <- Control_vs_HxP[order(Control_vs_HxP$pvalue),]
head(Control_vs_HxP)

# Heat vs HeatxPercipitation
Heat_vs_HxP <- get_results(dds,"Treatment","Heat","HeatxPercipitation")
summary(Heat_vs_HxP)
Heat_vs_HxP <- Heat_vs_HxP[order(Heat_vs_HxP$pvalue),]
head(Heat_vs_HxP)

# E+ vs E-
Endo_vs_No_Edno <- get_results(dds,"Endophyte","Negative","Positive")
summary(Endo_vs_No_Edno)
Endo_vs_No_Edno <- Endo_vs_No_Edno[order(Endo_vs_No_Edno$pvalue),]
head(Endo_vs_No_Edno)


#######################
# Volcano Plot Function 
#######################
# Define the function to create a volcano plot
create_volcano_plot <- function(results, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Volcano Plot") {
  
  # Ensure the results have the required columns
  if(!all(c("log2FoldChange", "padj") %in% colnames(results))) {
    stop("Results must contain 'log2FoldChange' and 'padj' columns.")
  }
  
  # Create a column to indicate significance
  results$significance <- "Not Significant"
  results$significance[results$padj < pvalue_threshold & results$log2FoldChange > log2FC_threshold] <- "Significant Upregulated"
  results$significance[results$padj < pvalue_threshold & results$log2FoldChange < (log2FC_threshold * -1)] <- "Significant Downregulated"
  
  # Create the volcano plot
  p <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significance), alpha = 0.9, size = 3) +
    scale_color_manual(values = c("Significant Upregulated" = "red", "Not Significant" = "grey", "Significant Downregulated" = "blue")) +
    geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
    labs(x = "Log2 Fold Change", y = "-log10 Adjusted P-value") +
    ggtitle(title) +
    theme_minimal() + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 23),  # Increased title size
      axis.title = element_text(size = 21),               # Increased axis title size
      axis.text = element_text(size = 19),                # Increased axis label size
      legend.text = element_text(size = 19),              # Increased legend text size
      legend.title = element_text(size = 21))              # Increased legend title size 
  return(p)
}
#######################
# Volcano Plot Function different colored 
#######################
# Define the function to create a volcano plot
create_volcano_plot_2_dataset <- function(results1, results2, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Volcano Plot") {
  
  # Ensure the results have the required columns
  if(!all(c("log2FoldChange", "padj") %in% colnames(results1))) {
    stop("Results must contain 'log2FoldChange' and 'padj' columns.")
  }
significance <- rep("Not Significant", nrow(results2))
significance[results2$padj < pvalue_threshold & results2$log2FoldChange > log2FC_threshold] <- "Significant Upregulated"
significance[results2$padj < pvalue_threshold & results2$log2FoldChange < -log2FC_threshold] <- "Significant Downregulated"
  

# Ensure that results1 and results2 have the same set of rownames (gene identifiers)
if (!all(rownames(results1) %in% rownames(results2))) {
  stop("Some gene identifiers in results1 are not found in results2.")
}
if (!all(rownames(results2) %in% rownames(results1))) {
  stop("Some gene identifiers in results2 are not found in results1.")
}

matched_indices <- match(rownames(results1), rownames(results2))

results1$significance <- significance[matched_indices]

not_significant <- subset(results1, significance == "Not Significant")
upregulated <- subset(results1, significance == "Significant Upregulated")
downregulated <- subset(results1, significance == "Significant Downregulated")

  # Create the volcano plot
  p <- ggplot() +
    geom_point(data = not_significant, aes(x = log2FoldChange, y = -log10(padj)), color = "grey", size = 1, alpha = 1) +
    geom_point(data = upregulated, aes(x = log2FoldChange, y = -log10(padj)), color = "red", size = 3, alpha = 1) +
    geom_point(data = downregulated, aes(x = log2FoldChange, y = -log10(padj)), color = "blue", size = 3, alpha = 1) +
    geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
    labs(x = "Log2 Fold Change", y = "-log10 Adjusted P-value") +
    ggtitle(title) +
    theme_minimal() + 
    theme(      plot.title = element_text(hjust = 0.5, size = 23),  # Increased title size
                axis.title = element_text(size = 21),               # Increased axis title size
                axis.text = element_text(size = 19),                # Increased axis label size
                legend.text = element_text(size = 19),              # Increased legend text size
                legend.title = element_text(size = 21))              # Increased legend title size )
  return(p)
}



###############################
# Volcano Plots
###############################

volcano_plot_CvH <- create_volcano_plot(Control_vs_Heat, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Control vs Heat")
print(volcano_plot_CvH)

volcano_plot_CvHxP <- create_volcano_plot(Control_vs_HxP, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Control vs HeatxPercipitation")
print(volcano_plot_CvHxP)

volcano_plot_HvHxP<- create_volcano_plot(Heat_vs_HxP, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Heat vs HeatxPercipitation")
print(volcano_plot_HvHxP)

volcano_plot_EvNE <- create_volcano_plot(Endo_vs_No_Edno, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endophyte Positive vs Endohpyte Negative")
print(volcano_plot_EvNE)

###################################
# Heat maps
###################################

########## Variation caused by Treatments on all data  #########################
# 1. Extract normalized counts
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# Calculate variance for each gene
gene_vars <- apply(vsd_mat, 1, var)

# Select top 50 most variable genes
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:50]
vsd_sub <- vsd_mat[top_genes, ]

# 4. Create annotation for columns
ann_col <- data.frame(Clone = metadata$Clone)
rownames(ann_col) <- metadata$SampleName

# 5. Plot with pheatmap
pheatmap(vsd_sub,
         annotation_col = ann_col,
         scale = "row",              # normalize rows for better visual
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Top 50 DEGs")













