# This is R script used to create better graphs that analyze gene expression directly than those in 2a
# Because I only use harvest time in later analysis I drop month and year here in favor of only HravestTime
# Auther: Darrian Talamantes

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
library(variancePartition)
library(BiocGenerics) 
library(lme4)        
library(gridExtra)



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
###############################
# Running DeSeq2 
###############################
#Makes DeSeq data set
dds <- DESeqDataSetFromMatrix(countData = Featurecount,
                              colData = metadata,
                              design= ~ Clone + HarvestTime + Endophyte + Treatment)

#Run DeSeq function, 
dds_og <- DESeq(dds)
dds <- dds_og


# filter for genes that have 10 occurrences in 1/4 the samples
keep <- rowSums(counts(dds) >= 5) >= (ncol(dds) / 4)
dds <- dds[keep, ]

####################################
# Function to get results.
####################################
get_results <- function(DeSeqData,column,t1,t2){
  if(!is.character(column) || !is.character(t1) || !is.character(t2)) {
    stop("The column and treatment levels (t1, t2) must be provided as strings.")
  }
  output <- results(DeSeqData, contrast = c(column, t1, t2))
  output <- output[order(output$padj),]
  return(output)
}



###################################
# Heat maps
###################################

########## Variation caused by Treatments on all data  #########################
# 1. Extract normalized counts
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# Calculate variance for each gene
gene_vars <- apply(vsd_mat, 1, var)

# Select top 50 most variable genes
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:50]
vsd_sub <- vsd_mat[top_genes, ]

# 4. Create annotation for columns
ann_col <- data.frame(Clone = metadata$Clone)
rownames(ann_col) <- metadata$SampleName

# 5. Plot with pheatmap
pheatmap(vsd_sub,
         annotation_col = ann_col,
         scale = "row",              # normalize rows for better visual
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Top 50 DEGs")


# Redo but with endophyte status
# 4. Create annotation for columns
ann_col <- data.frame(Endophyte = metadata$Endophyte)
rownames(ann_col) <- metadata$SampleName

# 5. Plot with pheatmap
pheatmap(vsd_sub,
         annotation_col = ann_col,
         scale = "row",              # normalize rows for better visual
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Top 50 DEGs")

# Redo one more time but with Time Point
# 4. Create annotation for columns
ann_col <- data.frame("Harvest Time" = metadata$HarvestTime)
rownames(ann_col) <- metadata$SampleName

# 5. Plot with pheatmap
pheatmap(vsd_sub,
         annotation_col = ann_col,
         scale = "row",              # normalize rows for better visual
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Top 50 DEGs")




# Basic PCA Plots
geno <- plotPCA(vsd, intgroup = "Clone")
HT <- plotPCA(vsd, intgroup = "HarvestTime")
treatment <- plotPCA(vsd, intgroup = "Treatment")
Endo <- plotPCA(vsd, intgroup = "Endophyte")

pca1 <- geno + labs(color = "Genotype") + theme_bw()
pca2 <- HT + labs(color = "Harvest Time") + theme_bw()
pca3 <- treatment + labs(color = "Treatment") + theme_bw()
pca4 <- Endo + labs(color = "Endophyte Status") + theme_bw()

grid.arrange(pca1, pca2, pca3, pca4, ncol = 2)

# Finding the variance of all varaibles 

meta <- as.data.frame(colData(dds))
form <- ~ (1|Clone) + (1|Treatment) + (1|HarvestTime) + (1|Endophyte) 
varPart <- fitExtractVarPartModel(vsd_mat, form, meta)
varPart_reordered <- varPart[, c("Clone", "HarvestTime", "Treatment", "Endophyte", "Residuals")]

p <- plotVarPart(varPart_reordered)
avg_var_explained <- colMeans(varPart_reordered)
print(avg_var_explained)
avg_df <- enframe(avg_var_explained, name = "Variable", value = "AvgVariance")

p + 
  geom_text(data = avg_df, aes(x = Variable, y = AvgVariance + 2, label = round(AvgVariance, 3)), 
            inherit.aes = FALSE, vjust = -35, size = 3.5,)


#  PCA Plots seperated by genotype


vsd_june <- vsd[, vsd$HarvestTime == "June_2016"]
vsd_june_CTE31 <- vsd_june[, vsd_june$Clone == "CTE31"]
vsd_june_CTE45 <- vsd_june[, vsd_june$Clone == "CTE45"]
vsd_june_CTE25 <- vsd_june[, vsd_june$Clone == "CTE25"]
vsd_june_CTE46 <- vsd_june[, vsd_june$Clone == "CTE46"]


junep1 <- plotPCA(vsd_june_CTE25, intgroup = "Treatment")
junep1_1 <- junep1 + labs(color = "Treatment") + theme_bw() + ggtitle("CTE25") +   theme(plot.title = element_text(hjust = 0.5))

junep2 <- plotPCA(vsd_june_CTE31, intgroup = "Treatment")
junep2_1 <- junep2 + labs(color = "Treatment") + theme_bw() + ggtitle("CTE31") +   theme(plot.title = element_text(hjust = 0.5))

junep3 <- plotPCA(vsd_june_CTE45, intgroup = "Treatment")
junep3_1 <- junep3 + labs(color = "Treatment") + theme_bw() + ggtitle("CTE45") +   theme(plot.title = element_text(hjust = 0.5))

june4 <- plotPCA(vsd_june_CTE46, intgroup = "Treatment")
junep4_1 <- june4 + labs(color = "Treatment") + theme_bw() + ggtitle("CTE46") +   theme(plot.title = element_text(hjust = 0.5))

junep1_1
junep2_1
junep3_1
junep4_1


small_stack <- grid.arrange(junep1_1, junep2_1, junep3_1, junep4_1, ncol = 2)
small_stack
grid.arrange(pca1, small_stack, nrow = 2, heights = c(1, 1.5))
