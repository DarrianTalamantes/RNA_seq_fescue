# Author: Darrian Talamantes
# This script is very similar to that of 2d
# but here when I filter the data to ensure that a gene appears at least 2 times 
# within any group. Ex: In total_degs_treatments I now keep any gene that has a count of
# 1 in any treatment group. Before I only kept genes with a count of 2. 
# This drastically changes the graphs. 



# NOTES: The filter on total_DEG data set and the filter on 
# _up and _down data sets are done differently! This may be problematic


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
library(cowplot)


# File locations
data_folder <- "/home/darrian/Documents/RNA_seq_fescue/r_data"
# Data recreation from previous plot
gene_counts = read.csv(paste0(data_folder, "/Epos_Eneg_Deseq2_contrast.csv"), header = TRUE) # from 2c


################################################################################
# Huge function that creates 4 data sets
# CTE all DEGs, CTE up or down regulated DEG count
# Treatments all DEGs, Treatments up or down regulated DEG count.
################################################################################

data_splitter <- function(GeneCount = gene_counts, cutoff = 2)
{
  # Prepare base data
  rownames(GeneCount) <- GeneCount$Gene
  data_only <- GeneCount[, -1]
  
  #### Subset per genotype first  ####
  cte_up <- data_only[, c("CTE25_Up", "CTE31_Up", "CTE45_Up", "CTE46_Up")]
  cte_down <- data_only[, c("CTE25_Down", "CTE31_Down", "CTE45_Down", "CTE46_Down")]
  
  # Now filter genes for CTE Up/Down separately
  cte_up_genes <- rownames(cte_up)[rowSums(cte_up >= cutoff) > 0]
  cte_down_genes <- rownames(cte_down)[rowSums(cte_down >= cutoff) > 0]
  cte_genes <- unique(c(cte_up_genes, cte_down_genes))
  
  # Subset just those genes
  cte_up_subset <- cte_up[cte_genes, , drop = FALSE]
  cte_down_subset <- cte_down[cte_genes, , drop = FALSE]
  
  # Convert rownames to Gene column for merging
  cte_up_subset$Gene <- rownames(cte_up_subset)
  cte_down_subset$Gene <- rownames(cte_down_subset)
  
  # Merge the two on Gene
  final_CTE_Up_Down <- merge(cte_up_subset, cte_down_subset, by = "Gene", all = TRUE)
  
  # Restore rownames and reorder columns
  rownames(final_CTE_Up_Down) <- final_CTE_Up_Down$Gene
  final_CTE_Up_Down$Gene <- NULL
  
  #### Repeat for Treatment  ####
  treat_up <- data_only[, c("Heat_Up", "Control_Up", "HeatxPercipitation_Up")]
  treat_down <- data_only[, c("Heat_Down", "Control_Down", "HeatxPercipitation_Down")]
  
  # Identify genes to retain
  treat_up_genes <- rownames(treat_up)[rowSums(treat_up >= cutoff) > 0]
  treat_down_genes <- rownames(treat_down)[rowSums(treat_down >= cutoff) > 0]
  treat_genes <- unique(c(treat_up_genes, treat_down_genes))
  
  # Subset the relevant genes
  treat_up_subset <- treat_up[treat_genes, , drop = FALSE]
  treat_down_subset <- treat_down[treat_genes, , drop = FALSE]
  
  # Add Gene column for merge
  treat_up_subset$Gene <- rownames(treat_up_subset)
  treat_down_subset$Gene <- rownames(treat_down_subset)
  
  # Merge by gene
  final_Treats_Up_Down <- merge(treat_up_subset, treat_down_subset, by = "Gene", all = TRUE)
  
  # Restore rownames and reorder
  rownames(final_Treats_Up_Down) <- final_Treats_Up_Down$Gene
  final_Treats_Up_Down$Gene <- NULL
  
  #### This will combine the up and down for count of all DEGs  ####
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
  
  total_degs_genos <- subset(final_DEGs, select = c ("CTE25", "CTE31", "CTE45", "CTE46"))
  total_degs_treatments <- subset(final_DEGs, select = c ("Heat", "Control", "HeatxPercipitation"))
  
  
  # Filteres out any genes that are not found to have a count of cutoff in at least 1 group
  top_genes_list2 <- apply(total_degs_genos, 2, function(col) {
    filtered_col <- col[col >= cutoff]
  })
  
  total_degs_genos$Gene <- rownames(total_degs_genos)
  # Combine all top gene names into one unique set
  top_genes_unique2 <- unique(unlist(lapply(top_genes_list2, names)))
  # Subset the original table with only those genes
  total_degs_genos <- total_degs_genos[total_degs_genos$Gene %in% top_genes_unique2, ]
  total_degs_genos$Gene <- NULL
  
  
  # Filteres out any genes that are not found to have a count of cutoff in at least 1 group
  top_genes_list3 <- apply(total_degs_treatments, 2, function(col) {
    filtered_col <- col[col >= cutoff]
  })
  
  total_degs_treatments$Gene <- rownames(total_degs_treatments)
  # Combine all top gene names into one unique set
  top_genes_unique3 <- unique(unlist(lapply(top_genes_list3, names)))
  # Subset the original table with only those genes
  total_degs_treatments <- total_degs_treatments[total_degs_treatments$Gene %in% top_genes_unique3, ]
  total_degs_treatments$Gene <- NULL
  
  return(list(
    final_CTE_Up_Down = final_CTE_Up_Down,
    final_Treats_Up_Down = final_Treats_Up_Down,
    total_degs_genos = total_degs_genos,
    total_degs_treatments = total_degs_treatments
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

binary_data <- function(dataframe){
  binary_df <- dataframe %>%
    rownames_to_column(var = "Gene") %>%
    mutate(across(-Gene, ~ ifelse(. != 0, 1, 0)))
  return(binary_df)
}


#######################
# Using the functions
#######################
all_counts <- data_splitter(gene_counts, 1)

upsetter(all_counts$final_CTE_Up_Down, "E+ and E- DEGs by Genotype")
upsetter(all_counts$final_Treats_Up_Down, "E+ and E- DEGs by Treatment")

p1 <- upsetter(all_counts$total_degs_genos, "E+ and E- DEGs by Genotype No DEG Count Filter")
p2 <- upsetter(all_counts$total_degs_treatments, "E+ and E- DEGs by Treatment No DEG Count Filter")
plot_grid(p1, p2, labels = c("A", "B"), ncol = 1, align = "v")

binary_data(all_counts$final_CTE_Up_Down)

