# Objective: This file will take the significant go term files and extract their gene-id and locations
# it will then pull the relevant fasta files so you can use BLAST on them.


import pandas as pd
import re

# --- Input file paths ---
go_terms_file = "../Goatools_data/Control_Down_results.txt"
gtf_file = "/scratch/drt83172/Wallace_lab/RNA_SEQ/transcriptome_final/Fescue_transcriptome.gtf"
genome = "/scratch/drt83172/Wallace_lab/RNA_SEQ/Genome/Tall_fescue/tall_fescue_pv1.1.fasta"
output_file = "../Goatools_data/gene_locations.tsv"

# Read GO terms file
go_df = pd.read_csv(go_terms_file, sep="\t")

# Extract all unique gene IDs from the 'study_items' column
gene_ids = set()
for item in go_df['study_items']:
    gene_ids.update(g.strip() for g in item.split(", "))

print(f"Loaded {len(gene_ids)} unique gene IDs from GO terms file")

# Read GTF file
gtf_df = pd.read_csv(
    transcriptome_file,
    sep="\t",
    comment="#",
    header=None,
    names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
)

# Function to extract gene_id from the attribute column
def extract_gene_id(attr):
    for field in attr.split(";"):
        if field.strip().startswith("gene_id"):
            return field.split('"')[1]
    return None

gtf_df["gene_id"] = gtf_df["attribute"].apply(extract_gene_id)

# Filter for matching gene IDs
gtf_filtered = gtf_df[gtf_df["gene_id"].isin(gene_ids)][["gene_id", "seqname", "start", "end"]]

# Show a preview to console (first 10 rows)
print("\nExample matches found:")
for idx, row in gtf_filtered.head(10).iterrows():
    print(f"{row['gene_id']} -> {row['seqname']}:{row['start']}-{row['end']}")

# Save to file
gtf_filtered.to_csv(output_file, sep="\t", index=False)
print(f"\nSaved {len(gtf_filtered)} gene locations to {output_file}")