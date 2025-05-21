# This is R script uses the will partition the data then run deseq on it many different ways
# It is a very exploritory script.
# Best seperation was from results_list_CloneXHTxTreat which breaks everything down up excpet endophyte which can contrast
# Then makes stacked bar plot and a preliminary upset plot.

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
# Making Volcano Plots
################################################################################

plot_volcano <- function(res, title = "Volcano Plot") {
  res_df <- as.data.frame(res)
  res_df$significant <- with(res_df, padj < 0.05 & abs(log2FoldChange) > 1.5)
  
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significant)) +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    ggtitle(title) +
    xlab("Log2 Fold Change, E+ vs E-") +
    ylab("-Log10 adjusted p-value")
}

for (key in names(results_list_CloneXHTxTreat)) {
  dds <- results_list_CloneXHTxTreat[[key]]
  
  if (!inherits(dds, "DESeqDataSet")) {
    message("Skipping ", key, ": Not a DESeqDataSet")
    next
  }
  res <- results(dds, contrast = c("Endophyte", "Positive", "Negative"))  # change levels as needed
  p <- plot_volcano(res, title = paste("Volcano:", key))
  if (!is.null(p)) print(p)
}


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



################################################################################
# Stacked Bar Plot
################################################################################


gene_counts = read.csv(paste0(data_folder, "/Epos_Eneg_Deseq2_contrast.csv"), header = TRUE)

rownames(gene_counts) <- gene_counts$Gene
data_only <- gene_counts[, -1]

# Get top 10 genes for each column
top_genes_list <- apply(data_only, 2, function(col) {
  names(sort(col, decreasing = TRUE))[1:10]
})

# Combine all top gene names into one unique set
top_genes_unique <- unique(unlist(top_genes_list))

# Subset the original table with only those genes
top_genes_df <- gene_counts[gene_counts$Gene %in% top_genes_unique, ]

# View the result
top_genes_df

# Subset the data
top_genes_df_CTEup <- subset_df <- subset(top_genes_df, select = c("CTE25_Up", "CTE31_Up", "CTE45_Up", "CTE46_Up" ))
top_genes_df_CTEdown <- subset_df <- subset(top_genes_df, select = c("CTE25_Down", "CTE31_Down", "CTE45_Down", "CTE46_Down" ))
top_genes_df_CTEdown <- top_genes_df_CTEdown * -1

merged_df <- merge(top_genes_df_CTEup, top_genes_df_CTEdown, 
                   by = "row.names", 
                   all = TRUE)

#### Long data then plot
# Step 1: Move rownames to a column
merged_df$Gene <- merged_df$Row.names
merged_df$Row.names <- NULL
# Step 2: Pivot to long format
long_df <- merged_df %>%
  pivot_longer(-Gene, names_to = "CloneCondition", values_to = "Count") %>%
  mutate(Count = as.numeric(Count)) %>%         # Ensure Count is numeric
  filter(!is.na(Count))                         # Remove NA values


is.numeric(long_df$Count)

# Step 3: Extract clone names from column names (e.g., "CTE25_Up" → "CTE25")
long_df <- long_df %>%
  mutate(Clone = sub("_.*", "", CloneCondition))
long_df$Gene <- as.factor(long_df$Gene)


# Step 4: Plot 
# Explanation: This plot max is -4 and 4. Each clone can get a max of 1.
# This shows top 10 deregulated genes in all Clones
ggplot(long_df, aes(x = Gene, y = Count, fill = Clone)) +
  geom_bar(stat = "identity", position = "stack", colour="black") +
  scale_fill_manual(values=c("darkgoldenrod1", "darkslategray1","darkolivegreen2", "darkorchid3")) +
  theme_bw() +
  ylab("Count (Upregulated = +, Downregulated = -)") +
  xlab("Gene") +
  ggtitle("Differentially Expressed Genes of E+ and E- Within a Genotype") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  guides(fill=guide_legend(title="Genotype"))


################################################################################
# Upset Plots
################################################################################

# Data recreation from previous plot
gene_counts = read.csv(paste0(data_folder, "/Epos_Eneg_Deseq2_contrast.csv"), header = TRUE)



rownames(gene_counts) <- gene_counts$Gene
data_only <- gene_counts[, -1]

# Get top 10 genes for each column
top_genes_list <- apply(data_only, 2, function(col) {
  names(sort(col, decreasing = TRUE))[1:10]
})

# Combine all top gene names into one unique set
top_genes_unique <- unique(unlist(top_genes_list))

# Subset the original table with only those genes
top_genes_df <- gene_counts[gene_counts$Gene %in% top_genes_unique, ]

# View the result
top_genes_df

# Subset the data
top_genes_df_CTEup <- subset_df <- subset(top_genes_df, select = c("CTE25_Up", "CTE31_Up", "CTE45_Up", "CTE46_Up" ))
top_genes_df_CTEdown <- subset_df <- subset(top_genes_df, select = c("CTE25_Down", "CTE31_Down", "CTE45_Down", "CTE46_Down" ))
top_genes_df_CTEdown <- top_genes_df_CTEdown * -1

merged_df <- merge(top_genes_df_CTEup, top_genes_df_CTEdown, 
                   by = "row.names", 
                   all = TRUE)

# Binary data for upset plot
binary_df <- merged_df %>%
  mutate(across(-Gene, ~ ifelse(. != 0, 1, 0)))

upset(binary_df, 
      intersect = colnames(binary_df)[-which(names(binary_df) == "Gene")],
      name = "Genes",
      base_annotations = list('Intersection size' = intersection_size()))



























