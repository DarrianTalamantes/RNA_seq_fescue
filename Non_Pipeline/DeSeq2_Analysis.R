# This is R script used to analyze the genes within the RNA seq data
# Auther: Darrian Talamantes

# Notes:
# Before starting this script I deleted the suffix within the feature counts table

# Install Bioconductor if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

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
                              design= ~ Month + Clone + Year + Treatment + Endophyte, tidy = TRUE)

#Run DeSeq function, 
dds_og <- DESeq(dds)
 dds <- dds_og

# filter for genes that have 10 occurrences in 1/4 the samples
keep <- rowSums(counts(dds) >= 5) >= (ncol(dds) / 4)
dds <- dds[keep, ]

####################################
# Getting many comparisons
####################################
# Function to get results.
get_results <- function(DeSeqData,column,t1,t2){
  if(!is.character(column) || !is.character(t1) || !is.character(t2)) {
    stop("The column and treatment levels (t1, t2) must be provided as strings.")
  }
  output <- results(DeSeqData, contrast = c(column, t1, t2))
  output <- output[order(output$padj),]
  return(output)
}


# Control vs Heat
Control_vs_Heat <- get_results(dds,"Treatment","Control","Heat")
summary(Control_vs_Heat)
Control_vs_Heat <- Control_vs_Heat[order(Control_vs_Heat$pvalue),]
head(Control_vs_Heat)

# Control vs HeatxPercipitation
Control_vs_HxP <- get_results(dds,"Treatment","Control","HeatxPercipitation")
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


###############################
# Scatter Plot
###############################


# Scatter Plot
# par(mfrow=c(2,3))
# plotCounts(dds, gene="gene.117427.469.0", intgroup="Treatment")
# plotCounts(dds, gene="gene.94951.99.0", intgroup="Treatment")
# plotCounts(dds, gene="gene.122487.2597.0", intgroup="Treatment")
# plotCounts(dds, gene="gene.17870.375.2", intgroup="Treatment")
# plotCounts(dds, gene="gene.54571.320.0", intgroup="Treatment")
# plotCounts(dds, gene="gene.1127.808.0", intgroup="Treatment")

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
    geom_point(aes(color = significance), alpha = 0.9) +
    scale_color_manual(values = c("Significant Upregulated" = "red", "Not Significant" = "grey", "Significant Downregulated" = "blue")) +
    geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
    labs(x = "Log2 Fold Change", y = "-log10 Adjusted P-value", title = title) +
    theme_minimal()
  
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


