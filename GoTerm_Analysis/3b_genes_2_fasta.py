# Objective: This file will take the significant go term files and extract their gene-id and locations
# it will then pull the relevant fasta files so you can use BLAST on them.


import pandas as pd
import re

# --- Input file paths ---
go_file = "../Goatools_data/Control_Down_results.txt"
gtf_file = "/scratch/drt83172/Wallace_lab/RNA_SEQ/transcriptome_final/Fescue_transcriptome.gtf"
genome = "/scratch/drt83172/Wallace_lab/RNA_SEQ/Genome/Tall_fescue/tall_fescue_pv1.1.fasta"
output_file = "../Goatools_data/gene_locations.tsv"

# --- 1. Read GO file and collect all unique gene IDs ---
go_df = pd.read_csv(go_file, sep="\t")
all_genes = set()

for genes in go_df["study_items"]:
    for g in re.split(r",\s*", genes):
        all_genes.add(g.strip())

# --- 2. Read GTF file ---
gtf_cols = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
gtf = pd.read_csv(gtf_file, sep="\t", comment="#", names=gtf_cols)

# Extract gene_id from attributes
gtf["gene_id"] = gtf["attribute"].str.extract(r'gene_id "([^"]+)"')

# --- 3. Filter for matching gene IDs ---
gtf_filtered = gtf[gtf["gene_id"].isin(all_genes)]

# --- 4. Keep only needed columns & unique entries ---
gtf_filtered = gtf_filtered[["gene_id", "chr", "start", "end"]].drop_duplicates()

# --- 5. Save to file ---
gtf_filtered.to_csv(output_file, sep="\t", index=False)

print(f"Saved {len(gtf_filtered)} gene locations to {output_file}")