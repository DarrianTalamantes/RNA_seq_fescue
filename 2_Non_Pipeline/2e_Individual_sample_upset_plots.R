# This will look at individual samples within treatment groups
# This will tell us if its just one sample going off the rails


library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(grid)
library(data.table)
library(pheatmap)










################################################################################
# Making table comparing diff expression of E+ an E- within Clones or treatments.
################################################################################

#### Creating lists of different sample groups to look at
names(results_list_CloneXHTxTreat)
# Get all the names
all_names <- names(results_list_CloneXHTxTreat)
# Extract treatment and group
treatments <- sapply(strsplit(all_names, "_"), function(x) tail(x, 1))
names_by_treatment <- split(all_names, treatments)
# Extract by clone and group
Clones <- sapply(strsplit(all_names, "_"), function(x) head(x, 1))
names_by_clones <- split(all_names, Clones)
#View by treatment
length(names_by_treatment$Heat)
length(names_by_clones$CTE25)
#### This will do a contrast between E+ and E- of each sample dds object
#### Then label genes upregulated or downregulated
# Step 1: Create a list to hold results labeled with up/down/false
label_results <- list()

for (key in names(results_list_CloneXHTxTreat)) {
  dds <- results_list_CloneXHTxTreat[[key]]
  
  tryCatch({
    res <- results(dds, contrast = c("Endophyte", "Positive", "Negative"))
    res_df <- as.data.frame(res)
    res_df$Gene <- rownames(res_df)
    
    res_df$Label <- "FALSE"
    res_df$Label[res_df$padj < 0.05 & res_df$log2FoldChange >= 1.5] <- "Upregulated"
    res_df$Label[res_df$padj < 0.05 & res_df$log2FoldChange <= -1.5] <- "Downregulated"
    
    # Store as named vector
    label_vec <- setNames(res_df$Label, res_df$Gene)
    label_results[[key]] <- label_vec
  }, error = function(e) {
    message("Error processing ", key, ": ", e$message)
  })
}

# Step 2: Get all unique gene names
all_genes <- unique(unlist(lapply(label_results, names)))

# Step 3: Initialize empty data frame
summary_df <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)

# Step 4: Fill in each column with labels
for (key in names(label_results)) {
  vec <- label_results[[key]]
  # Match to all_genes, fill with FALSE where not present
  summary_df[[key]] <- vec[summary_df$Gene]
  summary_df[[key]][is.na(summary_df[[key]])] <- "FALSE"
}


#### This seperates the data into the groups and counts up the up and down
#### regulated genes by group.
# Create an output data frame to store counts
gene_counts <- data.frame(Gene = summary_df$Gene, stringsAsFactors = FALSE)

# Define treatment groups
treatment_lists <- list(
  Heat = names_by_treatment$Heat,
  Control = names_by_treatment$Control,
  HeatxPercipitation = names_by_treatment$HeatxPercipitation,
  CTE25 = names_by_clones$CTE25,
  CTE31 = names_by_clones$CTE31,
  CTE45 = names_by_clones$CTE45,
  CTE46 = names_by_clones$CTE46
  # Add more if needed
)

# Loop through treatments
for (treatment in names(treatment_lists)) {
  samples <- treatment_lists[[treatment]]
  # Ensure samples exist in summary_df
  samples <- intersect(samples, colnames(summary_df))
  print(length(samples))
  if (length(samples) == 0) next
  
  # Subset summary_df to relevant columns
  subset_df <- summary_df[, samples, drop = FALSE]
  
  # Count Upregulated
  up_count <- apply(subset_df, 1, function(x) sum(x == "Upregulated"))
  down_count <- apply(subset_df, 1, function(x) sum(x == "Downregulated"))
  
  # Add to final count table
  gene_counts[[paste0(treatment, "_Up")]] <- up_count 
  gene_counts[[paste0(treatment, "_Down")]] <- down_count 
}

# View final table
head(gene_counts)

# Save the data
write.csv(gene_counts,paste0(data_folder, "/Epos_Eneg_Deseq2_contrast.csv"), row.names = FALSE)






