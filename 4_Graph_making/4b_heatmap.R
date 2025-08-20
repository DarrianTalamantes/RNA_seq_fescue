# Purpose: This script groups the data into the smallest possibel groups
# thus all groups have the same genotype, harvest time, and treatment. It then runs
# a contrast on those groups. When counting the DEGs this script uses two factors to group them
# all DEG counts are grouped by genotype and Treatment.

# Load Libraries

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
# Making table comparing diff expression of E+ an E- within Clones or treatments.
################################################################################

#### Creating lists of different sample groups to look at
names(results_list_CloneXHTxTreat)
# Get all the names
all_names <- names(results_list_CloneXHTxTreat)

# Putting all DeSeq2 objects into groups
group_keys <- sub("^(CTE\\d+)_.*_([^_]+)$", "\\1_\\2", all_names)
names_by_geno_treat <- split(all_names, group_keys)

#View by treatment
length(names_by_geno_treat$CTE25_Control)
length(names_by_geno_treat$CTE46_Heat)
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

# Loop through treatments
for (treatment in names(names_by_geno_treat)) {
  samples <- names_by_geno_treat[[treatment]]
  # Ensure samples exist in summary_df
  samples <- intersect(samples, colnames(summary_df))
  print(length(samples))
  if (length(samples) == 0) next
  
  # Subset summary_df to relevant columns
  subset_df <- summary_df[, samples, drop = FALSE]
  
  # Count Upregulated
  up_count <- apply(subset_df, 1, function(x) sum(x == "Upregulated"))
  down_count <- apply(subset_df, 1, function(x) sum(x == "Downregulated"))
  deg_count <- apply(subset_df, 1, function(x) sum(x == "Downregulated" | x == "Upregulated"))
  # Add to final count table
  gene_counts[[paste0(treatment, "_Up")]] <- up_count 
  gene_counts[[paste0(treatment, "_Down")]] <- down_count 
  gene_counts[[paste0(treatment, "_DEGs")]] <- deg_count 
  
}

# View final table
head(gene_counts)


################################################################################
# Making Heat Map
################################################################################

gene_counts_DEGs <- gene_counts[ , grep("DEGs$", colnames(gene_counts))]
rownames(gene_counts_DEGs) <- summary_df$Gene
gene_counts_DEGs <- gene_counts_DEGs[rowSums(gene_counts_DEGs) != 0, ]

gene_count_DEGs_mat <- as.matrix(gene_counts_DEGs)


# order by treatment
treatment <- sub("^.*?_(.*?)_DEGs$", "\\1", colnames(gene_counts_DEGs))
order_idx <- order(treatment)
gene_count_DEGs_treatorder <- gene_counts_DEGs[ , order_idx]
gene_count_DEGs_treatorder_mat <- as.matrix(gene_count_DEGs_treatorder)

# order by treatment
treatment <- sub("^.*?_(.*?)_DEGs$", "\\1", colnames(gene_counts_DEGs))
order_idx <- order(treatment)
gene_count_DEGs_treatorder <- gene_counts_DEGs[ , order_idx]
gene_count_DEGs_treatorder_mat <- as.matrix(gene_count_DEGs_treatorder)


# order by genotype
genotype <- sub("^(CTE\\d+)_.*$", "\\1", colnames(gene_counts_DEGs))
order_idx <- order(genotype)
gene_count_DEGs_genoorder <- gene_counts_DEGs[ , order_idx]
gene_count_DEGs_genoorder_mat <- as.matrix(gene_count_DEGs_genoorder)


pheatmap(gene_count_DEGs_mat,
         scale = "row",           # scale each gene across samples
         show_rownames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Heatmap of DEG Count Within Same Genotype and Treatment",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

pheatmap(gene_count_DEGs_treatorder_mat,
         scale = "row",           # scale each gene across samples
         show_rownames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Heatmap of DEG Count Within Same Genotype and Treatment (Treatment Order)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

pheatmap(gene_count_DEGs_genoorder_mat,
         scale = "row",           # scale each gene across samples
         show_rownames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Heatmap of DEG Count Within Same Genotype and Treatment (Genotype Order)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))



################################################################################
# Heat map of Raw expression data
################################################################################

# Step 1 create lists of up/down regulated gene names and their expression data
# List to store both labels and expression
expr_list <- list()
label_results <- list()

for (key in names(results_list_CloneXHTxTreat)) {
  dds <- results_list_CloneXHTxTreat[[key]]
  
  tryCatch({
    res <- results(dds, contrast = c("Endophyte", "Positive", "Negative"))
    res_df <- as.data.frame(res)
    res_df$Gene <- rownames(res_df)
    
    # Label genes
    res_df$Label <- "FALSE"
    res_df$Label[res_df$padj < 0.05 & res_df$log2FoldChange >= 1.5] <- "Upregulated"
    res_df$Label[res_df$padj < 0.05 & res_df$log2FoldChange <= -1.5] <- "Downregulated"
    
    # Keep only up/downregulated genes
    degs <- res_df$Gene[res_df$Label %in% c("Upregulated", "Downregulated")]
    
    # Store the normalized counts for these genes
    norm_counts <- assay(vst(dds))   # or rlog(dds) if you prefer
    expr_list[[key]] <- norm_counts[degs, , drop=FALSE]  # rows=genes, cols=samples
    
    # Store labels
    label_results[[key]] <- setNames(res_df$Label[res_df$Gene %in% degs], degs)
    
  }, error = function(e) {
    message("Error processing ", key, ": ", e$message)
  })
}

# Step 2 : Get all unique DE genes across all datasets
all_de_genes <- unique(unlist(lapply(label_results, names)))


# Step 3: convert each dataset into long format with sample info
expr_list <- lapply(expr_list, function(x) {
  if (!is.data.frame(x)) as.data.frame(x) else x
})

long_list <- lapply(names(expr_list), function(dataset) {
  df <- expr_list[[dataset]]
  df$Gene <- rownames(df)   # keep gene ID
  df_long <- df %>%
    pivot_longer(
      cols = -Gene,
      names_to = "Sample",
      values_to = "Expression"
    ) %>%
    mutate(Dataset = dataset)  # keep track of dataset
  df_long
})

# Step 4: combine all datasets into one big long data frame
big_df <- bind_rows(long_list)

# Step 5: optionally filter to DE genes only
big_df <- big_df %>% filter(Gene %in% all_de_genes)

big_df_filtered <- big_df %>% 
  filter(Gene %in% all_de_genes)

head(big_df_filtered)

# Step 6: Heatmap creation
heatmap_matrix <- big_df_filtered %>%
  select(Sample, Gene, Expression) %>%
  pivot_wider(names_from = Gene, values_from = Expression, values_fill = 0) %>%
  as.data.frame()

# Keep Sample names as rownames
rownames(heatmap_matrix) <- heatmap_matrix$Sample
heatmap_matrix$Sample <- NULL

# Create row annotation with Genotype
sample_annotation <- big_df_filtered %>%
  select(Sample, Dataset) %>%
  distinct() %>%
  mutate(Genotype = sub("_.*", "", Dataset)) %>%  # Extract text before first "_"
  column_to_rownames("Sample")

# Define colors for Genotypes
genotypes <- unique(sample_annotation$Genotype)
genotype_colors <- setNames(RColorBrewer::brewer.pal(min(length(genotypes), 12), "Set3"),
                            genotypes)
ann_colors <- list(Genotype = genotype_colors)

# Plot heatmap
# Create row annotation with Genotype
sample_annotation <- big_df_filtered %>%
  select(Sample, Dataset) %>%
  distinct() %>%
  mutate(Genotype = sub("_.*", "", Dataset)) %>%  # Extract text before first "_"
  column_to_rownames("Sample")

# Define colors for Genotypes
genotypes <- unique(sample_annotation$Genotype)
genotype_colors <- setNames(RColorBrewer::brewer.pal(min(length(genotypes), 12), "Set3"),
                            genotypes)
ann_colors <- list(Genotype = genotype_colors)

# Plot
p <- pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_row = sample_annotation["Genotype"],
  annotation_colors = ann_colors,
  clustering_method = "complete",
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_row = FALSE
)

# Wrap heatmap in gtable
heatmap_grob <- p$gtable

# Create text grobs
y_label <- textGrob("Samples", rot=90, gp=gpar(fontsize=14))
x_label <- textGrob("Genes", gp=gpar(fontsize=14))

# Arrange using layout_matrix
layout_mat <- rbind(
  c(NA, 1),  # Top row: NA on left, X label on right
  c(2, 3)    # Bottom row: Y label on left, heatmap on right
)

grid.arrange(
  x_label,    # index 1 in layout
  y_label,    # index 2
  heatmap_grob, # index 3
  layout_matrix = layout_mat,
  widths = c(1, 8), heights = c(1, 8)
)


