# Purpose: This script will partition the data into groups that differ only by Endophyte
# and harvest time. These groups are the same in genotype and treatment. I then run
# a contrast via endophyte + and -. I then make a heatmap from the result.

# Load Libraries

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
# Function that partitions data by Clone and treatment
################################################################################





# partitions the data by Clone and Treatment.
dds_by_ClonexTreat <- function(CountsData, Metadata = metadata, CloneName, Treat){
  
  # Subset metadata
  meta_clone <- Metadata[Metadata$Clone == CloneName, ]
  meta_clone_endo <- meta_clone[meta_clone$Treatment == Treat,]
  
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
  message(sum(keep), " genes passed the filtering step.")
  
  
  return(dds_clone)
}

################################################################################
# Using the function to create datasets seperated by Clone and Treatment
################################################################################


# 1. Get all unique combinations of Clone and Treatment
combos <- expand.grid(Clone = unique(metadata$Clone),
                      Treatment = unique(metadata$Treatment),
                      stringsAsFactors = FALSE)

# 2. Loop over each combination and run the function
results_list_ClonexTreat <- list()

for (clone_name in unique(metadata$Clone)) {
  for (treat in unique(metadata$Treatment)) {
    
    result_key <- paste0(clone_name,"_",treat)
    
    tryCatch({
      res <- dds_by_ClonexTreat(Featurecount, metadata, clone_name, treat )
      results_list_ClonexTreat[[result_key]] <- res
      message("Successfully processed: ", result_key)
    }, error = function(e) {
      message("Error in: ", result_key, " - ", e$message)
    })
  }
}

results(results_list_ClonexTreat$CTE25_Control)
results(results_list_ClonexTreat$CTE46_HeatxPercipitation)


#### Need to add the rest of code from 2b or 2c


################################################################################
# Making table comparing diff expression of E+ an E- within Genotypes or treatments.
################################################################################

#### Creating lists of different sample groups to look at
names(results_list_ClonexTreat)
# Get all the names
all_names <- names(results_list_ClonexTreat)
# Extract treatment and group
treatments <- sapply(strsplit(all_names, "_"), function(x) tail(x, 1))
names_by_treatment <- split(all_names, treatments)
# Extract by clone and group
Genotypes <- sapply(strsplit(all_names, "_"), function(x) head(x, 1))
names_by_genotype <- split(all_names, Genotypes)
#View by treatment
length(names_by_treatment$Heat)
length(names_by_genotype$CTE25)
#### This will do a contrast between E+ and E- of each sample dds object
#### Then label genes upregulated or downregulated
# Step 1: Create a list to hold results labeled with up/down/false
label_results <- list()

for (key in names(results_list_ClonexTreat)) {
  dds <- results_list_ClonexTreat[[key]]
  
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

################################################################################
# heatmap making
################################################################################

# Filter out genes never marked as Up or Down regulated
summary_df_filtered <- summary_df[
  apply(summary_df[-1], 1, function(x) any(x != "FALSE")),
]

# Copy and convert categorical labels to numeric
heatmap_matrix <- summary_df_filtered
rownames(heatmap_matrix) <- heatmap_matrix$Gene
heatmap_matrix$Gene <- NULL

# Map text to numeric values
heatmap_matrix[] <- lapply(heatmap_matrix, function(x) {
  ifelse(x == "Upregulated", 1, ifelse(x == "Downregulated", -1, 0))
})

# Convert to matrix
heatmap_matrix <- as.matrix(heatmap_matrix)

# Plot heatmap
pheatmap(
  heatmap_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  main = "DEG Heatmap: Up/Downregulated"
)

# This heat map looks at samples in groups that all have the same genotype and treatment
# Thus samples differ by edophyte and harvest time. It runs a contrast on all these samples by 
# endophyte. It then labels upregulated and downregulated +1 and -1 respectivley and neither 0.
#This is what is graphed. 



################################################################################
# Upset plot
################################################################################

# Convert to binary: TRUE if Up or Down
binary_matrix <- summary_df
binary_matrix[] <- lapply(binary_matrix, function(x) {
  ifelse(x %in% c("Upregulated", "Downregulated"), 1, 0)
})

# Convert Gene to rownames for use
rownames(binary_matrix) <- summary_df$Gene
binary_matrix$Gene <- NULL

# Convert to data frame for ComplexUpset
binary_df <- as.data.frame(binary_matrix)
binary_df$Gene <- rownames(binary_matrix)

# Melt into long form
long_df <- pivot_longer(binary_df, -Gene, names_to = "Group", values_to = "Present")

# Spread back into wide form, but binary (gene by group)
wide_binary <- long_df %>%
  filter(Present == 1) %>%
  select(-Present) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Group, values_from = value, values_fill = 0)

# Plot
ComplexUpset::upset(wide_binary, intersect = names(wide_binary)[-1], base_annotations=list('Intersection size'=intersection_size()))

























