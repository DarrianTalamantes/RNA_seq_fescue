# Purpose: This script will partition the data into its smallest groups only differeing by E-
# and E+. This will then allow me to compare genotype by treatment effect

# Load Libraries

library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(grid)
library(data.table)
library(pheatmap)


# Load Data

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
# Function that partitions data by Clone,HarvestTime, and treatment
################################################################################

# partitions the data by Clone Harvest Time, and Treatment.
dds_by_CloneXHTxTreat <- function(CountsData, Metadata, CloneName, HT, Treat){
  
  # Subset metadata
  meta_clone <- Metadata[Metadata$Clone == CloneName, ]
  meta_clone_endo <- meta_clone[meta_clone$HarvestTime == HT,]
  meta_clone_endo <- meta_clone_endo[meta_clone_endo$Treatment == Treat,]
  
  # Subset counts to only include samples that are in the meta_clone_endo dataset
  counts <- CountsData[, colnames(CountsData) %in% meta_clone_endo$SampleName]
  
  # Reorder columns of counts to match row order of metadata
  counts <- counts[, match(meta_clone_endo$SampleName, colnames(counts))]
  # Now create dds for that subset
  dds_clone <- DESeqDataSetFromMatrix(countData = counts,
                                      colData = meta_clone_endo,
                                      design = ~ Endophyte)  # or other factors
  dds_clone <- DESeq(dds_clone)
  
  # filter for genes that have 10 occurrences in 1/4 the samples
  keep <- rowSums(counts(dds_clone) >= 5) >= (ncol(dds_clone) / 4)
  dds_clone <- dds_clone[keep, ]
  
  return(dds_clone)
}


################################################################################
# Using the function to create datasets seperated by Clone, HT, and Treatment
################################################################################


# 1. Get all unique combinations of Clone and Endophyte
combos <- expand.grid(Clone = unique(metadata$Clone),
                      HarvestTime = unique(metadata$HarvestTime),
                      Treatment = unique(metadata$Treatment),
                      stringsAsFactors = FALSE)

# 2. Loop over each combination and run the function
results_list_CloneXHTxTreat <- list()

for (clone_name in unique(metadata$Clone)) {
  for (HT in unique(metadata$HarvestTime)) {
    for (treat in unique(metadata$Treatment)) {
      
      
      result_key <- paste0(clone_name, "_", HT, "_",treat)
      
      tryCatch({
        res <- dds_by_CloneXHTxTreat(Featurecount, metadata, clone_name, HT, treat )
        results_list_CloneXHTxTreat[[result_key]] <- res
        message("Successfully processed: ", result_key)
      }, error = function(e) {
        message("Error in: ", result_key, " - ", e$message)
      })
    }
  }
}

results(results_list_CloneXHTxTreat$CTE46_October_2016_HeatxPercipitation)
results(results_list_CloneXHTxTreat$CTE46_June_2017_HeatxPercipitation)


################################################################################
# Making 3 tabels counting DEGs in every genotype.
################################################################################
#### Creating lists of different sample groups to look at
names(results_list_CloneXHTxTreat)
# Get all the names
all_names <- names(results_list_CloneXHTxTreat)

# Split by "_"
split_names <- strsplit(all_names, "_")

# Extract genotypes and treatments
genotypes <- unique(sapply(split_names, function(x) x[1]))
treatments <- unique(sapply(split_names, function(x) x[length(x)]))

# Create combinations
combinations <- expand.grid(genotypes = genotypes, treatments = treatments)

# Combine into a string
combinations$Combo <- paste(combinations$genotypes, combinations$treatments, sep = "_")

names_by_geno_n_treat <- split(all_names, combinations$Combo)


#### This will do a contrast between E+ and E- of each sample dds object
#### Then label genes upregulated or downregulated
#### label_results will be a data set of each group with its genes labeled FLASE upregulated or downregulated
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

# Loop through treatments
for (combo in names(names_by_geno_n_treat)) {
  samples <- names_by_geno_n_treat[[combo]]
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
  gene_counts[[paste0(combo, "_Up")]] <- up_count 
  gene_counts[[paste0(combo, "_Down")]] <- down_count 
  gene_counts[[paste0(combo, "_DEG")]] <- down_count +  up_count 
  
}

# View final table
head(gene_counts)


################################################################################
# Heat map!
################################################################################

gene_counts_mat <- as.matrix(subset(gene_counts, select = grepl("_DEG$", names(gene_counts))))

#Find rows with only 0

zero_var_rows <- apply(gene_counts_mat, 1, function(x) sd(x) == 0)
gene_counts_mat_filtered <- gene_counts_mat[!zero_var_rows, ]

gene_counts_mat_ordered <- gene_counts_mat_filtered[, order(colnames(gene_counts_mat_filtered))]

# Plot heatmap ordered by genotpye
pheatmap(gene_counts_mat_ordered,
         scale = "none",         # Normalize genes across conditions
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         fontsize_row = 8,
         fontsize_col = 10,
         labels_row = NA,
         color = colorRampPalette(c("white", "red"))(100),
         main = "DEG Heatmap Ordered by Genotype")


############ Heat map ordered by Treatment

col_names <- colnames(gene_counts_mat_filtered)
treatment_part <- sapply(strsplit(col_names, "_"), function(x) x[2])
ordered_col_indices <- order(treatment_part)
gene_counts_mat_ordered_2 <- gene_counts_mat_filtered[, ordered_col_indices]


pheatmap(gene_counts_mat_ordered_2,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         labels_row = NA,
         fontsize_row = 8,
         fontsize_col = 10,
         color = colorRampPalette(c("white", "red"))(100),
         main = "DEG Heatmap Ordered by Treatment")

