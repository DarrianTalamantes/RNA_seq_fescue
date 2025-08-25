# The point of this code is to create a heat map on the base expression data.
# I will subset the data into only the DEGs found in the 48 groups that are only different by E+ and E-
# I will make Y labels that are sorted first by genotype, then harvest time, then treatment, then endophyte




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
library(gridExtra)


################################################################################
# load in all data 
################################################################################

# Loading in Meta Data
metadata_loc <-"/home/darrian/Documents/RNA_seq_fescue/r_data/metadata_match_to_FeatureCounts.txt"
meta_data <- read.table(metadata_loc)


# Loading in Feature Counts
data_folder <- "/home/darrian/Documents/RNA_seq_fescue/r_data"
Featurecount_loc <- paste0(data_folder, "/feature_counts_name_fixed.txt")
Featurecount <- read.table(Featurecount_loc, header = TRUE)
Featurecount <- Featurecount[, colnames(Featurecount) %in% meta_data$SampleName]

ncol(Featurecount)
nrow(meta_data)

#Load in the list of DEGs
all_de_genes <- readRDS("/home/darrian/Documents/RNA_seq_fescue/r_data/all_de_genes.rds")
dds <- readRDS("/home/darrian/Documents/RNA_seq_fescue/r_data/dds.rds")
label_results <- readRDS("/home/darrian/Documents/RNA_seq_fescue/r_data/label_results")

################################################################################
# Prepare data
################################################################################

# 1: Transform counts
vst_data <- vst(dds)

# 2: Extract expression matrix
mat <- assay(vst_data)
mat <- mat[rownames(mat) %in% all_de_genes, ]

# 2.5: Subset genes to top 300
# gene_vars <- apply(mat, 1, var)
# top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:300]
# mat_subset <- mat[top_genes, ]
# mat <- mat_subset

# 3: Ordering the data
metadata_ordered <- meta_data[match(colnames(mat), meta_data$SampleName), ]
if(any(is.na(metadata_ordered$SampleName))){
  stop("Some samples in the matrix do not have matching metadata!")
}

# 4: Build a sorting key (Genotype → HarvestTime → Treatment → Endophyte)
metadata_ordered <- metadata_ordered %>%
  arrange(Genotype, HarvestTime, Treatment, Endophyte)

# 5: Reorder the matrix columns according to this sort
mat <- mat[, metadata_ordered$SampleName]

# 6: Create annotation data frame for pheatmap
annotation_col <- metadata_ordered[, c("Genotype", "HarvestTime", "Treatment", "Endophyte")]
rownames(annotation_col) <- metadata_ordered$SampleName

# 7: Graph it
breaks <- c(
  seq(0, 5, length.out = 50),     # blue → yellow
  seq(5.01, 10, length.out = 50), # yellow → orange
  seq(10.01, 15, length.out = 50) # orange → red
)

# Generate the same number of colors - 1
my_colors <- colorRampPalette(c("blue", "yellow", "orange", "red"))(length(breaks) - 1)

# Plot
p <- pheatmap(
  mat,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "black",
  color = my_colors,
  breaks = breaks
)

# Wrap heatmap in gtable
heatmap_grob <- p$gtable

# Create text grobs
y_label <- textGrob("Samples", rot=90, gp=gpar(fontsize=14))
x_label <- textGrob("Genes", gp = gpar(fontsize = 14),hjust = 0.5, x = unit(0.4, "npc")  
)
# Arrange using layout_matrix
layout_mat <- rbind(
  c(NA, 1),  # Top row: NA on left, X label on right
  c(2, 3)    # Bottom row: Y label on left, heatmap on right
)

heat <- grid.arrange(
  x_label,    # index 1 in layout
  y_label,    # index 2
  heatmap_grob, # index 3
  layout_matrix = layout_mat,
  widths = c(1, 8), heights = c(1, 8)
)


# ggsave(heat,file ="4d_Normalized_count_heatmap.png", width = 10, height = 10, dpi = 300 )



################################################################################
# Heat map of normalized count data only DEG within a group
################################################################################


# Step 1: Find the group of each sample
sample_to_group <- sapply(colnames(mat), function(samp) {
  # Genotype
  geno <- str_extract(samp, "^CTE\\d+")
  
  # Harvest time
  harvest_code <- str_extract(samp, "O\\d+|J\\d+")
  harvest <- case_when(
    harvest_code == "O17" ~ "October_2017",
    harvest_code == "O16" ~ "October_2016",
    harvest_code == "J16" ~ "June_2016",
    harvest_code == "J17" ~ "June_2017",
    TRUE ~ harvest_code
  )
  
  # Treatment (C, H, HP) – fixed
  parts <- str_split(samp, "_")[[1]]
  treat_code <- str_extract(parts[3], "HP|H|C")
  treat <- case_when(
    treat_code == "HP" ~ "HeatxPercipitation",
    treat_code == "H"  ~ "Heat",
    treat_code == "C"  ~ "Control",
    TRUE ~ "Other"
  )
  
  paste(geno, harvest, treat, sep = "_")
})



# Make a copy of mat to modify
mat_filtered <- mat

# Loop over samples
for (samp in colnames(mat_filtered)) {
  
  group <- sample_to_group[samp]           # get the group for this sample
  deg_genes <- names(label_results[[group]])    # genes that are DEGs in this group
  
  # Set expression to 0 for genes not in deg_genes
  mat_filtered[!rownames(mat_filtered) %in% deg_genes, samp] <- 0
}


breaks <- c(
  seq(0, .1, length.out = 50),     # blue → yellow
  seq(.2, 5, length.out = 50), # yellow → orange
  seq(5.1, 10, length.out = 50) # orange → red
)

# Generate the same number of colors - 1
my_colors <- colorRampPalette(c("blue", "yellow", "orange", "red"))(length(breaks) - 1)


p <- pheatmap(
  mat_filtered,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "black",
  color = my_colors,
  breaks = breaks
)

# Wrap heatmap in gtable
heatmap_grob <- p$gtable

# Create text grobs
y_label <- textGrob("Samples", rot=90, gp=gpar(fontsize=14))
x_label <- textGrob("Genes", gp = gpar(fontsize = 14),hjust = 0.5, x = unit(0.4, "npc")  
)
# Arrange using layout_matrix
layout_mat <- rbind(
  c(NA, 1),  # Top row: NA on left, X label on right
  c(2, 3)    # Bottom row: Y label on left, heatmap on right
)

heat2 <- grid.arrange(
  x_label,    # index 1 in layout
  y_label,    # index 2
  heatmap_grob, # index 3
  layout_matrix = layout_mat,
  widths = c(1, 8), heights = c(1, 8)
)


# ggsave(heat2,file ="4d_Normalized_count_heatmap_2.png", width = 10, height = 10, dpi = 300 )

### Getting DEG count per group

# Convert to binary: anything >0 becomes 1
mat_binary <- (mat_filtered > 0) * 1  
col_sums <- colSums(mat_binary)

col_sums
mean(col_sums)



#### Checking Stuff (Can delete later) #########

all_zero_genes <- rownames(mat_filtered)[rowSums(mat_filtered != 0) == 0]

if (length(all_zero_genes) == 0) {
  cat("All genes have at least one non-zero value.\n")
} else {
  cat("Some genes are all zeros:\n")
  print(all_zero_genes)
}

num_nonzero_genes <- sum(rowSums(mat_filtered != 0) > 0)
num_nonzero_genes






