# Author: Darrian Talamantes
# This script will use the gene counts made in 2c from the partitioned data "results_list_CloneXHTxTreat" 
# with the contrast of E+ and E-. These counts are then used to make various different upset plots.
# split by treatment or genotype, then up and down regulated. or simply all DEGs counted up 


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

# Data recreation from previous plot
gene_counts = read.csv(paste0(data_folder, "/Epos_Eneg_Deseq2_contrast.csv"), header = TRUE) # from 2c
# File locations
data_folder <- "/home/darrian/Documents/RNA_seq_fescue/r_data"

################################################################################
# Huge function that creates 4 data sets
# CTE all DEGs, CTE up or down regulated DEG count
# Treatments all DEGs, Treatments up or down regulated DEG count.
################################################################################

data_splitter <- function(GeneCount = gene_counts, cutoff = 2)
{
  rownames(GeneCount) <- GeneCount$Gene
  data_only <- GeneCount[, -1]
  
  top_genes_list <- apply(data_only, 2, function(col) {
    # Filter out values that occur fewer than 2 times
    filtered_col <- col[col >= cutoff]
  })
  # Combine all top gene names into one unique set
  top_genes_unique <- unique(unlist(lapply(top_genes_list, names)))
  # Subset the original table with only those genes
  top_genes_df <- GeneCount[GeneCount$Gene %in% top_genes_unique, ]
  # View the result
  
  top_genes_df
  # This will count up all the Up and Down seperatly for CTE  
  top_genes_df_CTEup <- subset(top_genes_df, select = c("CTE25_Up", "CTE31_Up", "CTE45_Up", "CTE46_Up" ))
  top_genes_df_CTEdown <- subset(top_genes_df, select = c("CTE25_Down", "CTE31_Down", "CTE45_Down", "CTE46_Down" ))
  top_genes_df_CTEdown <- top_genes_df_CTEdown
  top_genes_df_CTEup$Gene <- rownames(top_genes_df_CTEup)
  top_genes_df_CTEdown$Gene <- rownames(top_genes_df_CTEdown)
  final_CTE_Up_Down <- merge(top_genes_df_CTEup, top_genes_df_CTEdown, 
                             by = "Gene", all = TRUE)
  rownames(final_CTE_Up_Down) <- final_CTE_Up_Down$Gene
  final_CTE_Up_Down$Gene <- NULL
  
  # This will count up all the Up and Down seperatly for Treatment  
  top_genes_df_Treatsup <- subset(top_genes_df, select = c("Heat_Up", "Control_Up", "HeatxPercipitation_Up" ))
  top_genes_df_Treatsdown <- subset(top_genes_df, select = c("Heat_Down", "Control_Down", "HeatxPercipitation_Down"))
  top_genes_df_Treatsdown <- top_genes_df_Treatsdown
  top_genes_df_Treatsup$Gene <- rownames(top_genes_df_Treatsup)
  top_genes_df_Treatsdown$Gene <- rownames(top_genes_df_Treatsdown)
  final_Treats_Up_Down <- merge(top_genes_df_Treatsup, top_genes_df_Treatsdown, 
                                by = "Gene", all = TRUE)
  rownames(final_Treats_Up_Down) <- final_Treats_Up_Down$Gene
  final_Treats_Up_Down$Gene <- NULL
  
  # This will combine the up and down for count of all DEGs  
  ids <- unique(sub("_.*", "", colnames(data_only)))
  sums <- list()
  for (cte in ids) {
    cols <- grep(paste0("^", cte, "_"), colnames(data_only), value = TRUE)
    message("Working on ", cte, " using columns: ", paste(cols, collapse = ", "))
    
    if (length(cols) > 0) {
      sums[[cte]] <- rowSums(data_only[, cols, drop = FALSE])
    } else {
      warning("No columns found for ", cte)
    }
  }
  
  final_DEGs <- as.data.frame(sums)
  rownames(final_DEGs) <- rownames(data_only)
  
  return(list(
    top_genes_df = top_genes_df,
    final_CTE_Up_Down = final_CTE_Up_Down,
    final_Treats_Up_Down = final_Treats_Up_Down,
    final_DEGs = final_DEGs
  ))
  
}



################################################################################
# Function to get data ready for upset plot and to make upset plot
################################################################################
# Binary data for upset plot

upsetter <- function(dataframe, title_name = "DEFAULT TITLE"){
  
  binary_df <- dataframe %>%
    rownames_to_column(var = "Gene") %>%
    mutate(across(-Gene, ~ ifelse(. != 0, 1, 0)))
  
  
  plot1 <- upset(binary_df,
                 intersect = colnames(binary_df)[-which(names(binary_df) == "Gene")],
                 name = "Unique Genes",
                 base_annotations = list('Intersection size' = intersection_size())) +
                 labs(title = title_name)

return(plot1)
}


#######################
# Using the functions
#######################
all_counts <- data_splitter(gene_counts, 2)
total_degs <- all_counts$final_DEGs
total_degs_genos <- subset(total_degs, select = c ("CTE25", "CTE31", "CTE45", "CTE46"))
total_degs_genos <- total_degs_genos[rowSums(total_degs_genos) != 0, ]
total_degs_treatments <- subset(total_degs, select = c ("Heat", "Control", "HeatxPercipitation"))
total_degs_treatments <- total_degs_treatments[rowSums(total_degs_treatments) != 0, ]


upsetter(all_counts$final_CTE_Up_Down, "E+ and E- DEGs by Genotype")
upsetter(all_counts$final_Treats_Up_Down, "E+ and E- DEGs by Treatment")

upsetter(total_degs_genos, "E+ and E- DEGs by Genotype")
upsetter(total_degs_treatments, "E+ and E- DEGs by Treatment")







