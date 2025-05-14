# This is R script partitions the data by clone and endophyte status then runs DeSeq2
# Author: Darrian Talamantes

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





################################################################################
# Function that partitions data by Endophtye and Clone
################################################################################

# partitions the data by CLone and Endophyte and returns deseq object
dds_by_CloneXEndo <- function(CountsData, Metadata, CloneName, EndoStatus){
  
  # Subset metadata for one Clone (e.g., Clone1)
  meta_clone <- Metadata[Metadata$Clone == CloneName, ]
  meta_clone_endo <- meta_clone[meta_clone$Endophyte == EndoStatus,]
  
  # Subset counts to only include samples that are in the meta_clone_endo dataset
  counts <- CountsData[, colnames(CountsData) %in% meta_clone_endo$SampleName]
  
  # Reorder columns of counts to match row order of metadata
  counts <- counts[, match(meta_clone_endo$SampleName, colnames(counts))]
  # Now create dds for that subset
  dds_clone <- DESeqDataSetFromMatrix(countData = counts,
                                       colData = meta_clone_endo,
                                       design = ~ Month + Year + Treatment)  # or other factors
  dds_clone <- DESeq(dds_clone)
  
  # filter for genes that have 10 occurrences in 1/4 the samples
  keep <- rowSums(counts(dds_clone) >= 5) >= (ncol(dds_clone) / 4)
  dds_clone <- dds_clone[keep, ]
  
  return(dds_clone)
}

################################################################################
# Function that partitions data by CLone only
################################################################################

# make deseq object that is based only on Clone name
dds_by_Clone <-  function(CountsData, Metadata, CloneName) {
  
  # Subset metadata for one Clone (e.g., Clone1)
  meta_clone <- Metadata[Metadata$Clone == CloneName, ]
  
  # Subset counts to only include samples that are in the meta_clone dataset
  counts <- CountsData[, colnames(CountsData) %in% meta_clone$SampleName]
  
  # Reorder columns of counts to match row order of metadata
  counts <- counts[, match(meta_clone$SampleName, colnames(counts))]
  
  # Create DESeq2 object
  dds_clone <- DESeqDataSetFromMatrix(countData = counts,
                                      colData = meta_clone,
                                      design = ~ Month + Year + Endophyte + Treatment)  # or include other factors if needed
  
  dds_clone <- DESeq(dds_clone)
  
  # Filter for genes with counts >=5 in at least 1/4 of samples
  keep <- rowSums(counts(dds_clone) >= 5) >= (ncol(dds_clone) / 4)
  dds_clone <- dds_clone[keep, ]
  
  return(dds_clone)
}


################################################################################
# Function that partitions data by Clone and Month
################################################################################

# partitions the data by CLone and Month and returns deseq object
dds_by_CloneXMonth <- function(CountsData, Metadata, CloneName, month){
  
  # Subset metadata for one Clone (e.g., Clone1)
  meta_clone <- Metadata[Metadata$Clone == CloneName, ]
  meta_clone_endo <- meta_clone[meta_clone$Month == month,]
  
  # Subset counts to only include samples that are in the meta_clone_endo dataset
  counts <- CountsData[, colnames(CountsData) %in% meta_clone_endo$SampleName]
  
  # Reorder columns of counts to match row order of metadata
  counts <- counts[, match(meta_clone_endo$SampleName, colnames(counts))]
  # Now create dds for that subset
  dds_clone <- DESeqDataSetFromMatrix(countData = counts,
                                      colData = meta_clone_endo,
                                      design = ~ Year + Endophyte + Treatment)  # or other factors
  dds_clone <- DESeq(dds_clone)
  
  # filter for genes that have 10 occurrences in 1/4 the samples
  keep <- rowSums(counts(dds_clone) >= 5) >= (ncol(dds_clone) / 4)
  dds_clone <- dds_clone[keep, ]
  
  return(dds_clone)
}

################################################################################
# Function that partitions data by Clone and HarvestTime
################################################################################

# partitions the data by CLone and Month and returns deseq object
dds_by_CloneXHT <- function(CountsData, Metadata, CloneName, HT){
  
  # Subset metadata for one Clone (e.g., Clone1)
  meta_clone <- Metadata[Metadata$Clone == CloneName, ]
  meta_clone_endo <- meta_clone[meta_clone$HarvestTime == HT,]
  
  # Subset counts to only include samples that are in the meta_clone_endo dataset
  counts <- CountsData[, colnames(CountsData) %in% meta_clone_endo$SampleName]
  
  # Reorder columns of counts to match row order of metadata
  counts <- counts[, match(meta_clone_endo$SampleName, colnames(counts))]
  # Now create dds for that subset
  dds_clone <- DESeqDataSetFromMatrix(countData = counts,
                                      colData = meta_clone_endo,
                                      design = ~ Endophyte + Treatment)  # or other factors
  dds_clone <- DESeq(dds_clone)
  
  # filter for genes that have 10 occurrences in 1/4 the samples
  keep <- rowSums(counts(dds_clone) >= 5) >= (ncol(dds_clone) / 4)
  dds_clone <- dds_clone[keep, ]
  
  return(dds_clone)
}

################################################################################
# Using the function to create datasets seperated by Clone and Endo
################################################################################


# 1. Get all unique combinations of Clone and Endophyte
combos <- expand.grid(Clone = unique(metadata$Clone),
                      Endophyte = unique(metadata$Endophyte),
                      stringsAsFactors = FALSE)

# 2. Loop over each combination and run the function
results_list_CloneXEndo <- list()

for (clone_name in unique(metadata$Clone)) {
  for (endo_status in unique(metadata$Endophyte)) {
    
    result_key <- paste0(clone_name, "_", endo_status)
    
    tryCatch({
      res <- dds_by_CloneXEndo(Featurecount, metadata, clone_name, endo_status)
      results_list_CloneXEndo[[result_key]] <- res
      message("Successfully processed: ", result_key)
    }, error = function(e) {
      message("Error in: ", result_key, " - ", e$message)
    })
  }
}
results(results_list_CloneXEndo$CTE25_Negative)
results(results_list_CloneXEndo$CTE46_Positive)


################################################################################
# Using the function to create datasets seperated by Clone and Month
################################################################################


# 1. Get all unique combinations of Clone and Month
combos <- expand.grid(Clone = unique(metadata$Clone),
                      Endophyte = unique(metadata$Month),
                      stringsAsFactors = FALSE)

# 2. Loop over each combination and run the function
results_list_by_ClonexMonth <- list()

for (clone_name in unique(metadata$Clone)) {
  for (month in unique(metadata$Month)) {
    
    result_key <- paste0(clone_name, "_", month)
    
    tryCatch({
      res <- dds_by_CloneXMonth(Featurecount, metadata, clone_name, month)
      results_list_by_ClonexMonth[[result_key]] <- res
      message("Successfully processed: ", result_key)
    }, error = function(e) {
      message("Error in: ", result_key, " - ", e$message)
    })
  }
}
results(results_list_by_ClonexMonth$CTE25_June)
results(results_list_by_ClonexMonth$CTE46_October)

################################################################################
# Using the function to create datasets seperated by Clone and HarvestTime
################################################################################


# 1. Get all unique combinations of Clone and Harvest Time
combos <- expand.grid(Clone = unique(metadata$Clone),
                      Endophyte = unique(metadata$HarvestTime),
                      stringsAsFactors = FALSE)

# 2. Loop over each combination and run the function
results_list_by_ClonexHT <- list()

for (clone_name in unique(metadata$Clone)) {
  for (HT in unique(metadata$HarvestTime)) {
    
    result_key <- paste0(clone_name, "_", HT)
    
    tryCatch({
      res <- dds_by_CloneXHT(Featurecount, metadata, clone_name, HT)
      results_list_by_ClonexHT[[result_key]] <- res
      message("Successfully processed: ", result_key)
    }, error = function(e) {
      message("Error in: ", result_key, " - ", e$message)
    })
  }
}
results(results_list_by_ClonexHT$CTE46_October_2017)
results(results_list_by_ClonexHT$CTE46_June_2016)

################################################################################
# Using the function to create datasets separated only by clone
################################################################################

# 1. Create a list to hold the results
results_list_by_clone <- list()

for (clone_name in unique(metadata$Clone)) {
  try({
    result <- dds_by_Clone(Featurecount, metadata, clone_name)
    results_list_by_clone[[clone_name]] <- result
  }, silent = TRUE)
}

results(results_list_by_clone$CTE25)
results(results_list_by_clone$CTE46)

################################################################################
# Create heatmap for data that is split by comparison and dendrogram seperaator
################################################################################

generate_heatmap_pheatmap <- function(comparison_name = "CTE25_Negative", results_list, metadata, seperator = "Treatment") {
  if (!is.character(comparison_name) || length(comparison_name) != 1) {
    stop("comparison_name must be a single string.")
  }
  dds <- results_list[[comparison_name]]
  
  # Extract normalized counts
  vsd <- vst(dds, blind = FALSE)
  vsd_mat <- assay(vsd)
  
  # Calculate variance for each gene
  gene_vars <- apply(vsd_mat, 1, var)
  
  # Select top 50 most variable genes
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:50]
  vsd_sub <- vsd_mat[top_genes, ]
  
  # Annotation for columns (samples)
  ann_col <- metadata %>%
    dplyr::filter(SampleName %in% colnames(vsd_sub)) %>%
    dplyr::select(SampleName, seperator)
  rownames(ann_col) <- ann_col$SampleName
  ann_col <- ann_col[seperator, drop = FALSE]
  
  # Plot
  pheatmap(vsd_sub,
           annotation_col = ann_col,
           scale = "row",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           show_rownames = TRUE,
           show_colnames = TRUE,
           main = paste("Top 50 Variable Genes:", comparison_name))
}

################################################################################
# generating final phenotypes
################################################################################
colnames(metadata)
names(results_list_CloneXEndo)
names(results_list_by_clone)
names(results_list_by_ClonexMonth)

generate_heatmap_pheatmap("CTE25_Negative",results_list,metadata,"Month")
generate_heatmap_pheatmap("CTE25_Positive",results_list,metadata,"Month")
generate_heatmap_pheatmap("CTE25",results_list_by_clone,metadata,"Month")

generate_heatmap_pheatmap("CTE45_Negative",results_list,metadata,"HarvestTime")
generate_heatmap_pheatmap("CTE45_Positive",results_list,metadata,"HarvestTime")
generate_heatmap_pheatmap("CTE45",results_list_by_clone,metadata,"HarvestTime")

generate_heatmap_pheatmap("CTE45_June",results_list_by_ClonexMonth,metadata,"Endophyte")
generate_heatmap_pheatmap("CTE45_October",results_list_by_ClonexMonth,metadata,"Treatment")
generate_heatmap_pheatmap("CTE25_June",results_list_by_ClonexMonth,metadata,"HarvestTime")
generate_heatmap_pheatmap("CTE25_October",results_list_by_ClonexMonth,metadata,"HarvestTime")





0colnames(metadata)


























































# 
# 
# 
# ################################################################################
# # DeSeq
# ################################################################################
# # remaking og dds file
# dds <- dds_og
# # Grouping our treatment groups
# dds$group <- factor(paste0(dds$Endophyte, dds$Treatment))
# design(dds) <- ~ group
# 
# #Run DeSeq function, 
# dds <- DESeq(dds)
# resultsNames(dds) # This will show you all your comparisions you want to look at
# 
# # filter for genes that have 10 occurrences in 1/4 the samples
# keep <- rowSums(counts(dds) >= 5) >= (ncol(dds) / 4)
# dds <- dds[keep, ]
# 
# ####################################
# # Getting contrats of interactions
# ####################################
# resultsNames(dds) # This will show you all your comparisions you want to look at
# 
# # The Treatments #
# 
# # Endo negative, heat v control
# EndoNeg_HeatxControl <- get_results(dds,"group","NegativeHeat","NegativeControl")
# summary(EndoNeg_HeatxControl)
# EndoNeg_HeatxControl <- EndoNeg_HeatxControl[order(EndoNeg_HeatxControl$pvalue),]
# head(EndoNeg_HeatxControl)
# 
# # Endo positive, heat v control
# EndoPos_HeatxControl <- get_results(dds,"group","PositiveHeat","PositiveControl")
# summary(EndoPos_HeatxControl)
# EndoPos_HeatxControl <- EndoPos_HeatxControl[order(EndoPos_HeatxControl$pvalue),]
# head(EndoPos_HeatxControl)
# 
# # Endo positive, HP v control
# EndoPos_HPxControl <- get_results(dds,"group","PositiveHeatxPercipitation","PositiveControl")
# summary(EndoPos_HPxControl)
# EndoPos_HPxControl <- EndoPos_HPxControl[order(EndoPos_HPxControl$pvalue),]
# head(EndoPos_HPxControl)
# 
# # Endo negative, HP v control
# EndoNeg_HPxControl <- get_results(dds,"group","NegativeHeatxPercipitation","NegativeControl")
# summary(EndoNeg_HPxControl)
# EndoNeg_HPxControl <- EndoNeg_HPxControl[order(EndoNeg_HPxControl$pvalue),]
# head(EndoNeg_HPxControl)
# 
# # Heat x Heat Percipitation #
# 
# # Endo positive, HP v control
# EndoPos_HPxHeat <- get_results(dds,"group","PositiveHeatxPercipitation","PositiveHeat")
# summary(EndoPos_HPxHeat)
# EndoPos_HPxHeat <- EndoPos_HPxHeat[order(EndoPos_HPxHeat$pvalue),]
# head(EndoPos_HPxHeat)
# 
# # Endo negative, HP v control
# EndoNeg_HPxHeat <- get_results(dds,"group","NegativeHeatxPercipitation","NegativeHeat")
# summary(EndoNeg_HPxHeat)
# EndoNeg_HPxHeat <- EndoNeg_HPxHeat[order(EndoNeg_HPxHeat$pvalue),]
# head(EndoNeg_HPxHeat)
# 
# 
# # Endophyte negative vs Positive #
# 
# # Endo Positive v negative , control
# Control_NegxPos <- get_results(dds,"group","NegativeControl","PositiveControl")
# summary(Control_NegxPos)
# Control_NegxPos <- Control_NegxPos[order(Control_NegxPos$pvalue),]
# head(Control_NegxPos)
# 
# # Endo Positive v negative , Heat
# Heat_NegxPos <- get_results(dds,"group","NegativeHeat","PositiveHeat")
# summary(Heat_NegxPos)
# Heat_NegxPos <- Heat_NegxPos[order(Heat_NegxPos$pvalue),]
# head(Heat_NegxPos)
# 
# # Endo Positive v negative , PxH
# PxH_NegxPos <- get_results(dds,"group","NegativeHeatxPercipitation","PositiveHeatxPercipitation")
# summary(PxH_NegxPos)
# PxH_NegxPos <- PxH_NegxPos[order(PxH_NegxPos$pvalue),]
# head(PxH_NegxPos)
# 
# ###################################
# # Volcano plots of interactions
# ###################################
# separator <- ggplot() + theme_void() + 
#   theme(panel.background = element_rect(fill = "black", colour = "black"))
# # Heat x Control
# volcano_plot_EndoNeg_HeatxControl <- create_volcano_plot(EndoNeg_HeatxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "E-,  Heat x Control")
# print(volcano_plot_EndoNeg_HeatxControl)
# 
# volcano_plot_EndoPos_HeatxControl <- create_volcano_plot(EndoPos_HeatxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "E+,  Heat x Control")
# print(volcano_plot_EndoPos_HeatxControl)
# 
# combined_plot1 <- ggarrange(volcano_plot_EndoNeg_HeatxControl, separator, volcano_plot_EndoPos_HeatxControl, ncol = 3, nrow =1, common.legend = TRUE, legend = "bottom", widths = c(1 ,.005 ,1))
# annotate_figure(combined_plot1, top = text_grob("Heat x Control Volcano Plots", face = "bold", size = 24))
# 
# # HeatxPresipitation x control
# volcano_plot_EndoPos_HPxControl <- create_volcano_plot(EndoPos_HPxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "E+, Heat with Percipitation x Control")
# print(volcano_plot_EndoPos_HPxControl)
# 
# volcano_plot_EndoNeg_HPxControl <- create_volcano_plot(EndoNeg_HPxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "E-, Heat with Percipitation x Control")
# print(volcano_plot_EndoNeg_HPxControl)
# 
# combined_plot2 <- ggarrange(volcano_plot_EndoNeg_HPxControl, separator, volcano_plot_EndoPos_HPxControl, ncol = 3, nrow =1, common.legend = TRUE, legend = "bottom", widths = c(1 ,.005 ,1))
# annotate_figure(combined_plot2, top = text_grob("Heat with Percipitation x Control Volcano Plots", face = "bold", size = 24))
# 
# # HeatxPresipitation x Heat
# volcano_plot_EndoPos_HPxHeat <- create_volcano_plot(EndoPos_HPxHeat, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "E+, Heat with Percipitation x Heat")
# print(volcano_plot_EndoPos_HPxHeat)
# 
# volcano_plot_EndoNeg_HPxHeat <- create_volcano_plot(EndoNeg_HPxHeat, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "E-, Heat with Percipitation x Heat")
# print(volcano_plot_EndoNeg_HPxHeat)
# 
# combined_plot3 <- ggarrange(volcano_plot_EndoNeg_HPxHeat, separator, volcano_plot_EndoPos_HPxHeat, ncol = 3, nrow =1, common.legend = TRUE, legend = "bottom", widths = c(1 ,.005 ,1))
# annotate_figure(combined_plot3, top = text_grob("Heat with Percipitation x Heat Volcano Plots", face = "bold", size = 24))
# 
# # Treatments, Negative x Positive
# volcano_plot_Heat_NegxPos <- create_volcano_plot(Heat_NegxPos, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Heat treatment, E- x E+")
# print(volcano_plot_Heat_NegxPos)
# 
# volcano_plot_HP_NegxPos <- create_volcano_plot(PxH_NegxPos, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Heat with Percipitation, E- x E+")
# print(volcano_plot_HP_NegxPos)
# 
# volcano_plot_Control_NegxPos <- create_volcano_plot(Control_NegxPos, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Control, Endo Negactive x Positive")
# print(volcano_plot_Control_NegxPos)
# 
# combined_plot4 <- ggarrange(volcano_plot_Heat_NegxPos, separator, volcano_plot_HP_NegxPos, ncol = 3, nrow =1, common.legend = TRUE, legend = "bottom", widths = c(1 ,.005 ,1))
# annotate_figure(combined_plot4, top = text_grob("Treatments by Endophyte Status", face = "bold", size = 24))
# 
# 
# 
# #######################################
# # PCA Plots function
# #######################################
# make_PCA <- function(dds, PCx = "PC1", PCy = "PC2", colors = "Treatment" ){
#   
#   vsd <- vst(dds, blind = FALSE)
#   # Doing a regular PCA without auto making it via DeSeq2
#   vst_data <- assay(vsd) # creates an expression matrix
#   pca_result <- prcomp(t(vst_data), center = TRUE, scale. = TRUE) # does PCA
#   explained_variance <- summary(pca_result)$importance[2, 1:100] * 100 # extracting explained variance
#   pca_data <- as.data.frame(pca_result$x) #savong as data frame
#   pca_data$sample <- rownames(pca_data)
#   sample_info <- as.data.frame(colData(dds))
#   pca_data <- merge(pca_data, sample_info, by.x = "sample", by.y = "row.names") # adding metadata
#   
#   PCx_sym <- sym(PCx)
#   PCy_sym <- sym(PCy)
#   colors_sym <- sym(colors)
#   
#   pca_plot <- ggplot(pca_data, aes(x = !!PCx_sym, y = !!PCy_sym, color = !!colors_sym)) + 
#     geom_point(size = 4) +
#     labs(
#       title = "PCA Plot",
#       x = paste(PCx, " - ", round(explained_variance[1], 1), "%", sep=""),
#       y = paste(PCy, " - ", round(explained_variance[2], 1), "%", sep="")
#     ) +
#     theme_minimal() +
#     theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
#   
#   return(pca_plot)
# }
# 
# #######################################
# # PCA Plots
# #######################################
# 
# vsd <- vst(dds, blind = FALSE)
# plotPCA(vsd, intgroup=c("Clone"))
# 
# # Doing a regular PCA without auto making it via DeSeq2
# vst_data <- assay(vsd) # creates an expression matrix
# pca_result <- prcomp(t(vst_data), center = TRUE, scale. = TRUE) # does PCA
# explained_variance <- summary(pca_result)$importance[2, 1:100] * 100 # extracting explained variance
# pca_data <- as.data.frame(pca_result$x) #savong as data frame
# pca_data$sample <- rownames(pca_data)
# sample_info <- as.data.frame(colData(dds))
# pca_data <- merge(pca_data, sample_info, by.x = "sample", by.y = "row.names") # adding metadata
# 
# pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Time )) + 
#   geom_point(size = 4) +
#   labs(
#     title = "PCA Plot",
#     x = paste("PC1 - ", round(explained_variance[1], 1), "%", sep=""),
#     y = paste("PC2 - ", round(explained_variance[2], 1), "%", sep="")
#   ) +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
# print(pca_plot)
# 
# # PC 2 seems to be Month
# # Pc 3 seems to be year
# 
# #### I think we can investigate more. We have to remove CTE25 and CTE31. This will allow us to better see differences since those are missing months in the data
# #### This analysis data with CTE25 and 31 removed ####
# # making filtered dds data
# MetaData_filtered <- subset(MetaData, (Clone %in% c("CTE25", "CTE31")))
# names_to_remove <- MetaData_filtered$Sample
# MetaData_filtered <- subset(MetaData, !(Clone %in% c("CTE25", "CTE31")))
# # Remove columns in df2 that match the names in df1
# Featurecount_filtered <- Featurecount[, !(colnames(Featurecount) %in% names_to_remove)]
# 
# dds_filtered <- DESeqDataSetFromMatrix(countData = Featurecount_filtered,
#                                        colData = MetaData_filtered,
#                                        design= ~ Month + Clone + Year + Endophyte + Treatment, tidy = TRUE)
# # filter for genes that have 10 occurrences in 1/4 the samples
# keep <- rowSums(counts(dds_filtered) >= 5) >= (ncol(dds_filtered) / 4)
# dds_filtered <- dds_filtered[keep, ]
# #PCA Make 
# # Proves Cone is PC1
# PCA_filtered <- make_PCA(dds_filtered,"PC1","PC2","Clone")
# print(PCA_filtered)
# # Proves PC2 and PC3 are the year and Month
# PCA_filtered <- make_PCA(dds_filtered,"PC2","PC3","Time")
# print(PCA_filtered)
# 
# 
# ########################################
# # Investigating the Clone
# ########################################
# 
# # # remaking og dds file
# # dds_clone <- dds_og
# # 
# # #Run DeSeq function, 
# # dds_clone <- DESeq(dds_clone)
# # 
# # # filter for genes that have 10 occurrences in 1/4 the samples
# # keep <- rowSums(counts(dds_clone) >= 5) >= (ncol(dds_clone) / 4)
# # dds_clone <- dds_clone[keep, ]
# # results(dds_clone)
# # # Clone 
# # CTE45xCTE46 <- get_results(dds,"Clone","CTE45","CTE46")
# # summary(CTE45xCTE46)
# # CTE45xCTE46 <- CTE45xCTE46[order(CTE45xCTE46$pvalue),]
# # head(CTE45xCTE46)
# # 
# # CTE45xCTE46 <- create_volcano_plot(CTE45xCTE46, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "CTE45 vs CTE46")
# # print(CTE45xCTE46)
# 
# ########################################
# # Creating list of significant genes
# ########################################
# 
# # Please run lines 267 - 328 before this to make the data
# ####### Function to get a significance table #########
# find_significance <- function(results, significant_table, log2FC_threshold = 2, pvalue_threshold = 0.05) {
#   
#   # Ensure the results have the required columns
#   if(!all(c("log2FoldChange", "padj") %in% colnames(results))) {
#     stop("Results must contain 'log2FoldChange' and 'padj' columns.")
#   }
#   
#   # getting the name of the dataframe as a string
#   name <- deparse(substitute(results))
#   
#   results$significance <- "Not Significant"
#   results$significance[results$padj < pvalue_threshold & results$log2FoldChange > log2FC_threshold] <- "Significant Upregulated"
#   results$significance[results$padj < pvalue_threshold & results$log2FoldChange < (log2FC_threshold * -1)] <- "Significant Downregulated"
#   
#   
#   df_subset <- results["significance"]
#   df_subset$Gene <- rownames(df_subset)
#   df_subset <- as.data.frame(df_subset)
#   merged_df <- merge(significant_table, df_subset, by = "Gene", all = TRUE)
#   # Here we would save the name of the datatable we are using as a string to be used in the next command.
#   colnames(merged_df)[colnames(merged_df) == "significance"] <- name
#   
#   return(merged_df)
# }
# 
# 
# ########## preparing to run significance function ###########
# #Make a list of datasets we want to analyze.  
# comparison_list <- c("EndoNeg_HeatxControl", "EndoPos_HeatxControl", "EndoPos_HPxControl", "EndoNeg_HPxControl", "EndoPos_HPxHeat", "EndoNeg_HPxHeat", "Control_NegxPos", "Heat_NegxPos", "PxH_NegxPos")
# gene_names <- rownames(dds)
# 
# # Creating blank data table to fill in with values
# significant_table <- data.table(matrix(0, nrow = length(gene_names), ncol = 1))
# significant_table <- as.data.frame(significant_table)
# rownames(significant_table) <- gene_names
# significant_table$Gene <- rownames(significant_table)
# head(significant_table)
# 
# ######## Running significance function ###########
# for (data_set in comparison_list) {
#   print(data_set)
#   data_frame <- get(data_set)
#   significant_table2 <- find_significance(data_frame, significant_table, 2, 0.05)
#   significant_table <- significant_table2
#   head(significant_table)
# }
# colnames(significant_table)[3:ncol(significant_table)] <- comparison_list
# head(significant_table)
# write.csv(significant_table, file = "significant_table.csv", row.names = FALSE)
# 
# 
# ################################################################################
# 
# ################################################################################
# EndoNeg_HeatxControl2 <- EndoNeg_HeatxControl
# 
# 
# EndoNeg_HeatxControl2$significance <- "Not Significant"
# EndoNeg_HeatxControl2$significance[EndoNeg_HeatxControl2$padj < .05 & EndoNeg_HeatxControl2$log2FoldChange > 2] <- "Significant Upregulated"
# EndoNeg_HeatxControl2$significance[EndoNeg_HeatxControl2$padj < .05 & EndoNeg_HeatxControl2$log2FoldChange < (2 * -1)] <- "Significant Downregulated"
# 
# 
# df_subset <- EndoNeg_HeatxControl2["significance"]
# df_subset$Gene <- rownames(df_subset)
# df_subset <- as.data.frame(df_subset)
# merged_df <- merge(significant_table, df_subset, by = "Gene", all = TRUE)
# # Here we would save the name of the datatable we are using as a string to be used in the next command.
# colnames(merged_df)[colnames(merged_df) == "significance"] <- "name1"
# 
# 
# significant_table2 <- find_significance(EndoNeg_HeatxControl, significant_table, 2, 0.05)
# significant_table <- significant_table2
# 
# significant_table2 <- find_significance(EndoPos_HeatxControl, significant_table, 2, 0.05)
# significant_table <- significant_table2
# 
# 
# head(significant_table)
# 
# class(data_name)
# head(data_name)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ############### Comparing 2 datasets using volcano plots #################
# 
# annotate_figure(combined_plot1, top = text_grob("Heat x Control Volcano Plots", face = "bold", size = 19))
# 
# create_volcano_plot_2_dataset(EndoNeg_HeatxControl, Heat_NegxPos, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Negative HeatxControl data colored by NegativexPositive")
# create_volcano_plot_2_dataset(Heat_NegxPos, EndoNeg_HeatxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Negative HeatxControl data colored by NegativexPositive")
# 
# create_volcano_plot_2_dataset(EndoPos_HeatxControl, Heat_NegxPos, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Positive HeatxControl data colored by NegativexPositive")
# create_volcano_plot_2_dataset(Heat_NegxPos, EndoPos_HeatxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Positive HeatxControl data colored by NegativexPositive")
# 
# 
# annotate_figure(combined_plot2, top = text_grob("Heat with Percipitation x Control Volcano Plots", face = "bold", size = 19))
# 
# create_volcano_plot_2_dataset(EndoNeg_HPxControl, PxH_NegxPos, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Negative HPxControl data colored by NegativexPositive")
# create_volcano_plot_2_dataset(PxH_NegxPos, EndoNeg_HPxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Negative HPxControl data colored by NegativexPositive")
# 
# create_volcano_plot_2_dataset(EndoPos_HPxControl, PxH_NegxPos, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Positive HPxControl data colored by NegativexPositive")
# create_volcano_plot_2_dataset(PxH_NegxPos, EndoPos_HPxControl, log2FC_threshold = 2, pvalue_threshold = 0.05, title = "Endo Positive HPxControl data colored by NegativexPositive")
# 
# 
# 
# 
