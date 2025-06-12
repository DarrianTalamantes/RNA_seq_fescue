# Notes: This script is supposed to run DEseq so I can compare heat stress and control 
# treatments. This will hopefully lead to lots more charachterized proteins.
# I am doing this because loooking at all the E+ and E- contrasts lead to many uncahrachterized proteins.

#Authoer: Darrian Talamantes 


#Loading Libraries

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
# Function that partitions data by Clone,HarvestTime, and Endophyte
################################################################################

# partitions the data by Clone Harvest Time, and Endophyte.
dds_by_CloneXHTxEndo <- function(CountsData, Metadata, CloneName, HT, Endo){
  
  # Subset metadata
  meta_clone <- Metadata[Metadata$Clone == CloneName, ]
  meta_clone_endo <- meta_clone[meta_clone$HarvestTime == HT,]
  meta_clone_endo <- meta_clone_endo[meta_clone_endo$Endophyte == Endo,]
  
  # Subset counts to only include samples that are in the meta_clone_endo dataset
  counts <- CountsData[, colnames(CountsData) %in% meta_clone_endo$SampleName]
  
  # Reorder columns of counts to match row order of metadata
  counts <- counts[, match(meta_clone_endo$SampleName, colnames(counts))]
  # Now create dds for that subset
  dds_clone <- DESeqDataSetFromMatrix(countData = counts,
                                      colData = meta_clone_endo,
                                      design = ~ Treatment)  # or other factors
  dds_clone <- DESeq(dds_clone)
  
  # filter for genes that have 10 occurrences in 1/4 the samples
  keep <- rowSums(counts(dds_clone) >= 5) >= (ncol(dds_clone) / 4)
  dds_clone <- dds_clone[keep, ]
  
  return(dds_clone)
}



################################################################################
# Using the function to create datasets seperated by Clone, HT, and Endophyte
################################################################################


# 1. Get all unique combinations
combos <- expand.grid(Clone = unique(metadata$Clone),
                      HarvestTime = unique(metadata$HarvestTime),
                      Endophyte = unique(metadata$Endophyte),
                      stringsAsFactors = FALSE)

# 2. Loop over each combination and run the function
results_list_CloneXHTxEndo <- list()

for (clone_name in unique(metadata$Clone)) {
  for (HT in unique(metadata$HarvestTime)) {
    for (treat in unique(metadata$Endophyte)) {
      
      
      result_key <- paste0(clone_name, "_", HT, "_",treat)
      
      tryCatch({
        res <- dds_by_CloneXHTxEndo(Featurecount, metadata, clone_name, HT, treat )
        results_list_CloneXHTxEndo[[result_key]] <- res
        message("Successfully processed: ", result_key)
      }, error = function(e) {
        message("Error in: ", result_key, " - ", e$message)
      })
    }
  }
}

results(results_list_CloneXHTxEndo$CTE25_June_2016_Negative)
results(results_list_CloneXHTxEndo$CTE25_June_2016_Positive)


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
# generating heatmaps
################################################################################

generate_heatmap_pheatmap("CTE25_June_2016_Positive",results_list_CloneXHTxEndo,metadata,"Treatment")
generate_heatmap_pheatmap("CTE25_June_2017_Positive",results_list_CloneXHTxEndo,metadata,"Treatment")





################################################################################
# Making table comparing diff expression of Heat vs Control within Clones or Endophyte status
################################################################################

#### Creating lists of different sample groups to look at
names(results_list_CloneXHTxEndo)
# Get all the names
all_names <- names(results_list_CloneXHTxEndo)
# Extract treatment and group
Endophytes <- sapply(strsplit(all_names, "_"), function(x) tail(x, 1))
names_by_Endophytes <- split(all_names, Endophytes)
# Extract by clone and group
Clones <- sapply(strsplit(all_names, "_"), function(x) head(x, 1))
names_by_clones <- split(all_names, Clones)
#View by treatment
length(names_by_Endophytes$Heat)
length(names_by_clones$CTE25)
#### This will do a contrast between E+ and E- of each sample dds object
#### Then label genes upregulated or downregulated
# Step 1: Create a list to hold results labeled with up/down/false
label_results <- list()

for (key in names(results_list_CloneXHTxEndo)) {
  dds <- results_list_CloneXHTxEndo[[key]]
  
  tryCatch({
    res <- results(dds, contrast = c("Treatment", "Heat", "Control"))
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
  Positive = names_by_Endophytes$Positive,
  Negative = names_by_Endophytes$Negative,
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
write.csv(gene_counts,paste0(data_folder, "/Heat_Control_Deseq2_contrast.csv"), row.names = FALSE)



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
  
  #### Repeat for Endophyte Status  ####
  enod_up <- data_only[, c("Positive_Up", "Negative_Up")]
  enod_down <- data_only[, c("Positive_Down", "Negative_Down")]
  
  # Identify genes to retain
  enod_up_genes <- rownames(enod_up)[rowSums(enod_up >= cutoff) > 0]
  enod_down_genes <- rownames(enod_down)[rowSums(enod_down >= cutoff) > 0]
  enod_genes <- unique(c(enod_up_genes, enod_down_genes))
  
  # Subset the relevant genes
  enod_up_subset <- enod_up[enod_genes, , drop = FALSE]
  enod_down_subset <- enod_down[enod_genes, , drop = FALSE]
  
  # Add Gene column for merge
  enod_up_subset$Gene <- rownames(enod_up_subset)
  enod_down_subset$Gene <- rownames(enod_down_subset)
  
  # Merge by gene
  final_enods_Up_Down <- merge(enod_up_subset, enod_down_subset, by = "Gene", all = TRUE)
  
  # Restore rownames and reorder
  rownames(final_enods_Up_Down) <- final_enods_Up_Down$Gene
  final_enods_Up_Down$Gene <- NULL
  
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
  total_degs_endophytes <- subset(final_DEGs, select = c ("Positive", "Negative"))
  
  
  # Filteres out any genes that are not found to have a count of 2 in at least 1 group
  top_genes_list2 <- apply(total_degs_genos, 2, function(col) {
    filtered_col <- col[col >= cutoff]
  })
  
  total_degs_genos$Gene <- rownames(total_degs_genos)
  # Combine all top gene names into one unique set
  top_genes_unique2 <- unique(unlist(lapply(top_genes_list2, names)))
  # Subset the original table with only those genes
  total_degs_genos <- total_degs_genos[total_degs_genos$Gene %in% top_genes_unique2, ]
  total_degs_genos$Gene <- NULL
  
  
  # Filteres out any genes that are not found to have a count of 2 in at least 1 group
  top_genes_list3 <- apply(total_degs_endophytes, 2, function(col) {
    filtered_col <- col[col >= cutoff]
  })
  
  total_degs_endophytes$Gene <- rownames(total_degs_endophytes)
  # Combine all top gene names into one unique set
  top_genes_unique3 <- unique(unlist(lapply(top_genes_list3, names)))
  # Subset the original table with only those genes
  total_degs_endophytes <- total_degs_endophytes[total_degs_endophytes$Gene %in% top_genes_unique3, ]
  total_degs_endophytes$Gene <- NULL
  
  return(list(
    final_CTE_Up_Down = final_CTE_Up_Down,
    final_enods_Up_Down = final_enods_Up_Down,
    total_degs_genos = total_degs_genos,
    total_degs_endophytes = total_degs_endophytes
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
all_counts <- data_splitter(gene_counts, 2)

upsetter(all_counts$final_CTE_Up_Down, "Heat vs Control DEGs by Genotype")
upsetter(all_counts$final_enods_Up_Down, "Heat vs Control DEGs by Endophyte Status")


upsetter(all_counts$total_degs_genos, "Heat vs Control DEGs by Genotype")
upsetter(all_counts$total_degs_endophytes, "Heat vs Control DEGs by Endophyte Status")


write.csv(all_counts$final_CTE_Up_Down,paste0(data_folder, "/Genotypes_Up_Down_reg_HeatvsControl.csv"), row.names = TRUE)
write.csv(all_counts$final_enods_Up_Down,paste0(data_folder, "/Endophytes_Up_Down_reg_HeatvsControl.csv"), row.names = TRUE)
