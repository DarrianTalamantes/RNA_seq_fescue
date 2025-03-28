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
library(unix)

# Increase memory of R to 12 GB
rlimit_as(1e12)

# File locations
MetaData_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/RNA_seq_fescue/Non_Pipeline/Meta_Data.csv"
Featurecount_loc <- "/home/darrian/Desktop/UGA/Wallace_Lab/Plant_Info/RNA_SEQ_backup_Data/feature_counts_namefix.txt"

# MetaData_loc <- "/home/drt06/Documents/Tall_fescue/RNA_seq_fescue/Non_Pipeline/Meta_Data.csv"
# Featurecount_loc <- "/home/drt06/Documents/Tall_fescue/Plant_Info/RNA_SEQ_backup_Data/feature_counts_namefix.txt"

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
                              design= ~ Month + Clone + Year + Endophyte + Treatment, tidy = TRUE)

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
# Heat map (Doesnt seem to work now)
###################################
# Check what comparison the data is doing. Seems like id have to rerun the data everytime for heatmaps. 
# I probably 
results(dds)

vsd <- vst(dds, blind = FALSE)
vsd_data <- assay(vsd)
# Select highly expressed genes
select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:100]
top200_data <- vsd_data[select, ]
# Convert the matrix to a long format for ggplot2
top200_long <- as.data.frame(top200_data) %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression")
# Add sample information
sample_info <- as.data.frame(colData(dds))
top200_long <- top200_long %>%
  left_join(sample_info, by = c("sample" = "Sample"))

# Create the heatmap ()
heatmap <- ggplot(Heat_vs_HxP, aes(x = sample, y = gene, fill = log2FoldChange)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),  # Remove x-axis labels
        axis.text.y = element_blank(),  # Remove y-axis labels
        panel.grid = element_blank()) +  # Remove grid lines
  labs(x = "Sample", y = "Gene", title = "Heatmap of Top 100 Gene Expression") +
  facet_grid(~ Treatment, scales = "free_x", space = "free_x") 
print(heatmap)




################################################################################
# DeSeq2 function interaction method
################################################################################
# remaking og dds file
dds <- dds_og
# Grouping our treatment groups
dds$group <- factor(paste0(dds$Endophyte, dds$Treatment))
design(dds) <- ~ group

#Run DeSeq function, 
dds <- DESeq(dds)
resultsNames(dds) # This will show you all your comparisions you want to look at

# filter for genes that have 10 occurrences in 1/4 the samples
keep <- rowSums(counts(dds) >= 5) >= (ncol(dds) / 4)
dds <- dds[keep, ]

####################################
# Getting contrats of interactions
####################################
resultsNames(dds) # This will show you all your comparisions you want to look at

# The Treatments #

# Endo negative, heat v control
EndoNeg_HeatxControl <- get_results(dds,"group","NegativeHeat","NegativeControl")
summary(EndoNeg_HeatxControl)
EndoNeg_HeatxControl <- EndoNeg_HeatxControl[order(EndoNeg_HeatxControl$pvalue),]
head(EndoNeg_HeatxControl)

# Endo positive, heat v control
EndoPos_HeatxControl <- get_results(dds,"group","PositiveHeat","PositiveControl")
summary(EndoPos_HeatxControl)
EndoPos_HeatxControl <- EndoPos_HeatxControl[order(EndoPos_HeatxControl$pvalue),]
head(EndoPos_HeatxControl)

# Endo positive, HP v control
EndoPos_HPxControl <- get_results(dds,"group","PositiveHeatxPercipitation","PositiveControl")
summary(EndoPos_HPxControl)
EndoPos_HPxControl <- EndoPos_HPxControl[order(EndoPos_HPxControl$pvalue),]
head(EndoPos_HPxControl)

# Endo negative, HP v control
EndoNeg_HPxControl <- get_results(dds,"group","NegativeHeatxPercipitation","NegativeControl")
summary(EndoNeg_HPxControl)
EndoNeg_HPxControl <- EndoNeg_HPxControl[order(EndoNeg_HPxControl$pvalue),]
head(EndoNeg_HPxControl)

# Heat x Heat Percipitation #

# Endo positive, HP v control
EndoPos_HPxHeat <- get_results(dds,"group","PositiveHeatxPercipitation","PositiveHeat")
summary(EndoPos_HPxHeat)
EndoPos_HPxHeat <- EndoPos_HPxHeat[order(EndoPos_HPxHeat$pvalue),]
head(EndoPos_HPxHeat)

# Endo negative, HP v control
EndoNeg_HPxHeat <- get_results(dds,"group","NegativeHeatxPercipitation","NegativeHeat")
summary(EndoNeg_HPxHeat)
EndoNeg_HPxHeat <- EndoNeg_HPxHeat[order(EndoNeg_HPxHeat$pvalue),]
head(EndoNeg_HPxHeat)


# Endophyte negative vs Positive #

# Endo Positive v negative , control
Control_NegxPos <- get_results(dds,"group","NegativeControl","PositiveControl")
summary(Control_NegxPos)
Control_NegxPos <- Control_NegxPos[order(Control_NegxPos$pvalue),]
head(Control_NegxPos)

# Endo Positive v negative , Heat
Heat_NegxPos <- get_results(dds,"group","NegativeHeat","PositiveHeat")
summary(Heat_NegxPos)
Heat_NegxPos <- Heat_NegxPos[order(Heat_NegxPos$pvalue),]
head(Heat_NegxPos)

# Endo Positive v negative , PxH
PxH_NegxPos <- get_results(dds,"group","NegativeHeatxPercipitation","PositiveHeatxPercipitation")
summary(PxH_NegxPos)
PxH_NegxPos <- PxH_NegxPos[order(PxH_NegxPos$pvalue),]
head(PxH_NegxPos)

###################################
# Volcano plots of interactions
###################################
separator <- ggplot() + theme_void() + 
  theme(panel.background = element_rect(fill = "black", colour = "black"))
# Heat x Control
volcano_plot_EndoNeg_HeatxControl <- create_volcano_plot(EndoNeg_HeatxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "E-,  Heat x Control")
print(volcano_plot_EndoNeg_HeatxControl)

volcano_plot_EndoPos_HeatxControl <- create_volcano_plot(EndoPos_HeatxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "E+,  Heat x Control")
print(volcano_plot_EndoPos_HeatxControl)

combined_plot1 <- ggarrange(volcano_plot_EndoNeg_HeatxControl, separator, volcano_plot_EndoPos_HeatxControl, ncol = 3, nrow =1, common.legend = TRUE, legend = "bottom", widths = c(1 ,.005 ,1))
annotate_figure(combined_plot1, top = text_grob("Heat x Control Volcano Plots", face = "bold", size = 24))

# HeatxPresipitation x control
volcano_plot_EndoPos_HPxControl <- create_volcano_plot(EndoPos_HPxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "E+, Heat with Percipitation x Control")
print(volcano_plot_EndoPos_HPxControl)

volcano_plot_EndoNeg_HPxControl <- create_volcano_plot(EndoNeg_HPxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "E-, Heat with Percipitation x Control")
print(volcano_plot_EndoNeg_HPxControl)

combined_plot2 <- ggarrange(volcano_plot_EndoNeg_HPxControl, separator, volcano_plot_EndoPos_HPxControl, ncol = 3, nrow =1, common.legend = TRUE, legend = "bottom", widths = c(1 ,.005 ,1))
annotate_figure(combined_plot2, top = text_grob("Heat with Percipitation x Control Volcano Plots", face = "bold", size = 24))

# HeatxPresipitation x Heat
volcano_plot_EndoPos_HPxHeat <- create_volcano_plot(EndoPos_HPxHeat, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "E+, Heat with Percipitation x Heat")
print(volcano_plot_EndoPos_HPxHeat)

volcano_plot_EndoNeg_HPxHeat <- create_volcano_plot(EndoNeg_HPxHeat, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "E-, Heat with Percipitation x Heat")
print(volcano_plot_EndoNeg_HPxHeat)

combined_plot3 <- ggarrange(volcano_plot_EndoNeg_HPxHeat, separator, volcano_plot_EndoPos_HPxHeat, ncol = 3, nrow =1, common.legend = TRUE, legend = "bottom", widths = c(1 ,.005 ,1))
annotate_figure(combined_plot3, top = text_grob("Heat with Percipitation x Heat Volcano Plots", face = "bold", size = 24))

# Treatments, Negative x Positive
volcano_plot_Heat_NegxPos <- create_volcano_plot(Heat_NegxPos, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Heat treatment, E- x E+")
print(volcano_plot_Heat_NegxPos)

volcano_plot_HP_NegxPos <- create_volcano_plot(PxH_NegxPos, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Heat with Percipitation, E- x E+")
print(volcano_plot_HP_NegxPos)

volcano_plot_Control_NegxPos <- create_volcano_plot(Control_NegxPos, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Control, Endo Negactive x Positive")
print(volcano_plot_Control_NegxPos)

combined_plot4 <- ggarrange(volcano_plot_Heat_NegxPos, separator, volcano_plot_HP_NegxPos, ncol = 3, nrow =1, common.legend = TRUE, legend = "bottom", widths = c(1 ,.005 ,1))
annotate_figure(combined_plot4, top = text_grob("Treatments by Endophyte Status", face = "bold", size = 24))



#######################################
# PCA Plots function
#######################################
make_PCA <- function(dds, PCx = "PC1", PCy = "PC2", colors = "Treatment" ){

vsd <- vst(dds, blind = FALSE)
# Doing a regular PCA without auto making it via DeSeq2
vst_data <- assay(vsd) # creates an expression matrix
pca_result <- prcomp(t(vst_data), center = TRUE, scale. = TRUE) # does PCA
explained_variance <- summary(pca_result)$importance[2, 1:100] * 100 # extracting explained variance
pca_data <- as.data.frame(pca_result$x) #savong as data frame
pca_data$sample <- rownames(pca_data)
sample_info <- as.data.frame(colData(dds))
pca_data <- merge(pca_data, sample_info, by.x = "sample", by.y = "row.names") # adding metadata

PCx_sym <- sym(PCx)
PCy_sym <- sym(PCy)
colors_sym <- sym(colors)

pca_plot <- ggplot(pca_data, aes(x = !!PCx_sym, y = !!PCy_sym, color = !!colors_sym)) + 
  geom_point(size = 4) +
  labs(
    title = "PCA Plot",
    x = paste(PCx, " - ", round(explained_variance[1], 1), "%", sep=""),
    y = paste(PCy, " - ", round(explained_variance[2], 1), "%", sep="")
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

return(pca_plot)
}

#######################################
# PCA Plots
#######################################

vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup=c("Clone"))

# Doing a regular PCA without auto making it via DeSeq2
vst_data <- assay(vsd) # creates an expression matrix
pca_result <- prcomp(t(vst_data), center = TRUE, scale. = TRUE) # does PCA
explained_variance <- summary(pca_result)$importance[2, 1:100] * 100 # extracting explained variance
pca_data <- as.data.frame(pca_result$x) #savong as data frame
pca_data$sample <- rownames(pca_data)
sample_info <- as.data.frame(colData(dds))
pca_data <- merge(pca_data, sample_info, by.x = "sample", by.y = "row.names") # adding metadata

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Time )) + 
  geom_point(size = 4) +
  labs(
    title = "PCA Plot",
    x = paste("PC1 - ", round(explained_variance[1], 1), "%", sep=""),
    y = paste("PC2 - ", round(explained_variance[2], 1), "%", sep="")
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
print(pca_plot)

# PC 2 seems to be Month
# Pc 3 seems to be year

#### I think we can investigate more. We have to remove CTE25 and CTE31. This will allow us to better see differences since those are missing months in the data
#### This analysis data with CTE25 and 31 removed ####
# making filtered dds data
MetaData_filtered <- subset(MetaData, (Clone %in% c("CTE25", "CTE31")))
names_to_remove <- MetaData_filtered$Sample
MetaData_filtered <- subset(MetaData, !(Clone %in% c("CTE25", "CTE31")))
# Remove columns in df2 that match the names in df1
Featurecount_filtered <- Featurecount[, !(colnames(Featurecount) %in% names_to_remove)]

dds_filtered <- DESeqDataSetFromMatrix(countData = Featurecount_filtered,
                              colData = MetaData_filtered,
                              design= ~ Month + Clone + Year + Endophyte + Treatment, tidy = TRUE)
# filter for genes that have 10 occurrences in 1/4 the samples
keep <- rowSums(counts(dds_filtered) >= 5) >= (ncol(dds_filtered) / 4)
dds_filtered <- dds_filtered[keep, ]
#PCA Make 
# Proves Cone is PC1
PCA_filtered <- make_PCA(dds_filtered,"PC1","PC2","Clone")
print(PCA_filtered)
# Proves PC2 and PC3 are the year and Month
PCA_filtered <- make_PCA(dds_filtered,"PC2","PC3","Time")
print(PCA_filtered)


########################################
# Investigating the Clone
########################################

# # remaking og dds file
# dds_clone <- dds_og
# 
# #Run DeSeq function, 
# dds_clone <- DESeq(dds_clone)
# 
# # filter for genes that have 10 occurrences in 1/4 the samples
# keep <- rowSums(counts(dds_clone) >= 5) >= (ncol(dds_clone) / 4)
# dds_clone <- dds_clone[keep, ]
# results(dds_clone)
# # Clone 
# CTE45xCTE46 <- get_results(dds,"Clone","CTE45","CTE46")
# summary(CTE45xCTE46)
# CTE45xCTE46 <- CTE45xCTE46[order(CTE45xCTE46$pvalue),]
# head(CTE45xCTE46)
# 
# CTE45xCTE46 <- create_volcano_plot(CTE45xCTE46, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "CTE45 vs CTE46")
# print(CTE45xCTE46)

########################################
# Creating list of significant genes
########################################

# Please run lines 267 - 328 before this to make the data
####### Function to get a significance table #########
find_significance <- function(results, significant_table, log2FC_threshold = 2, pvalue_threshold = 0.05) {
 
   # Ensure the results have the required columns
  if(!all(c("log2FoldChange", "padj") %in% colnames(results))) {
    stop("Results must contain 'log2FoldChange' and 'padj' columns.")
  }
  
  # getting the name of the dataframe as a string
  name <- deparse(substitute(results))
  
  results$significance <- "Not Significant"
  results$significance[results$padj < pvalue_threshold & results$log2FoldChange > log2FC_threshold] <- "Significant Upregulated"
  results$significance[results$padj < pvalue_threshold & results$log2FoldChange < (log2FC_threshold * -1)] <- "Significant Downregulated"
  
  
  df_subset <- results["significance"]
  df_subset$Gene <- rownames(df_subset)
  df_subset <- as.data.frame(df_subset)
  merged_df <- merge(significant_table, df_subset, by = "Gene", all = TRUE)
  # Here we would save the name of the datatable we are using as a string to be used in the next command.
  colnames(merged_df)[colnames(merged_df) == "significance"] <- name
  
  return(merged_df)
  }


########## preparing to run significance function ###########
#Make a list of datasets we want to analyze.  
comparison_list <- c("EndoNeg_HeatxControl", "EndoPos_HeatxControl", "EndoPos_HPxControl", "EndoNeg_HPxControl", "EndoPos_HPxHeat", "EndoNeg_HPxHeat", "Control_NegxPos", "Heat_NegxPos", "PxH_NegxPos")
gene_names <- rownames(dds)

# Creating blank data table to fill in with values
significant_table <- data.table(matrix(0, nrow = length(gene_names), ncol = 1))
significant_table <- as.data.frame(significant_table)
rownames(significant_table) <- gene_names
significant_table$Gene <- rownames(significant_table)
head(significant_table)

######## Running significance function ###########
for (data_set in comparison_list) {
  print(data_set)
  data_frame <- get(data_set)
  significant_table2 <- find_significance(data_frame, significant_table, 2, 0.05)
  significant_table <- significant_table2
  head(significant_table)
}
colnames(significant_table)[3:ncol(significant_table)] <- comparison_list
head(significant_table)
write.csv(significant_table, file = "significant_table.csv", row.names = FALSE)


################################################################################

################################################################################
EndoNeg_HeatxControl2 <- EndoNeg_HeatxControl


EndoNeg_HeatxControl2$significance <- "Not Significant"
EndoNeg_HeatxControl2$significance[EndoNeg_HeatxControl2$padj < .05 & EndoNeg_HeatxControl2$log2FoldChange > 2] <- "Significant Upregulated"
EndoNeg_HeatxControl2$significance[EndoNeg_HeatxControl2$padj < .05 & EndoNeg_HeatxControl2$log2FoldChange < (2 * -1)] <- "Significant Downregulated"


df_subset <- EndoNeg_HeatxControl2["significance"]
df_subset$Gene <- rownames(df_subset)
df_subset <- as.data.frame(df_subset)
merged_df <- merge(significant_table, df_subset, by = "Gene", all = TRUE)
# Here we would save the name of the datatable we are using as a string to be used in the next command.
colnames(merged_df)[colnames(merged_df) == "significance"] <- "name1"


significant_table2 <- find_significance(EndoNeg_HeatxControl, significant_table, 2, 0.05)
significant_table <- significant_table2

significant_table2 <- find_significance(EndoPos_HeatxControl, significant_table, 2, 0.05)
significant_table <- significant_table2


head(significant_table)

class(data_name)
head(data_name)
























############### Comparing 2 datasets using volcano plots #################

annotate_figure(combined_plot1, top = text_grob("Heat x Control Volcano Plots", face = "bold", size = 19))

create_volcano_plot_2_dataset(EndoNeg_HeatxControl, Heat_NegxPos, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Negative HeatxControl data colored by NegativexPositive")
create_volcano_plot_2_dataset(Heat_NegxPos, EndoNeg_HeatxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Negative HeatxControl data colored by NegativexPositive")

create_volcano_plot_2_dataset(EndoPos_HeatxControl, Heat_NegxPos, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Positive HeatxControl data colored by NegativexPositive")
create_volcano_plot_2_dataset(Heat_NegxPos, EndoPos_HeatxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Positive HeatxControl data colored by NegativexPositive")


annotate_figure(combined_plot2, top = text_grob("Heat with Percipitation x Control Volcano Plots", face = "bold", size = 19))

create_volcano_plot_2_dataset(EndoNeg_HPxControl, PxH_NegxPos, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Negative HPxControl data colored by NegativexPositive")
create_volcano_plot_2_dataset(PxH_NegxPos, EndoNeg_HPxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Negative HPxControl data colored by NegativexPositive")

create_volcano_plot_2_dataset(EndoPos_HPxControl, PxH_NegxPos, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Positive HPxControl data colored by NegativexPositive")
create_volcano_plot_2_dataset(PxH_NegxPos, EndoPos_HPxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Positive HPxControl data colored by NegativexPositive")

  
  

