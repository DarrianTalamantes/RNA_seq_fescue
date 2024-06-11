# This script uses R featurecounts for RNA seq pipline anlaysis.

args <- commandArgs(trailingOnly=TRUE)
bam_files <- strsplit(args[3], " ") # string to charachter vector
gtf_file <- args[2]
output_file <- args[1]

library(Rsubread)

# Run featureCounts
fc <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  nthreads = 8,
  isPairedEnd = TRUE
)

# Write counts to a file
write.table(fc$counts, file=output_file, sep="\t", quote=FALSE, col.names=NA)