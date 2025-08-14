# Objective: This file will take the significant go term files and extract their gene-id and locations
# it will then pull the relevant fasta files so you can use BLAST on them.
# Results: This program finds that out of the results files that mattered we can not find anything in our transcriptome file that matches.
# Remember Scallop2 creates the transcriptome with the gene.XX.XX.XX names
# 
 
import pandas as pd
import re

# --- Input file paths ---
go_terms_file = "../Goatools_data/Control_Down_results.txt"
transcriptome_file  = "/scratch/drt83172/Wallace_lab/RNA_SEQ/transcriptome_final/Fescue_transcriptome.gtf"
genome = "/scratch/drt83172/Wallace_lab/RNA_SEQ/Genome/Tall_fescue/tall_fescue_pv1.1.fasta"
output_file = "../Goatools_data/gene_locations.tsv"

# Read GO terms file
go_df = pd.read_csv(go_terms_file, sep="\t")

# Extract all unique gene IDs from the 'study_items' column
gene_ids = set()
for item in go_df['study_items']:
    gene_ids.update(g.strip() for g in item.split(", "))

print(f"Loaded {len(gene_ids)} unique gene IDs from GO terms file")

# Show a few extracted gene IDs for sanity check
print("\nExample gene IDs from GO terms file:")
for g in list(gene_ids)[:10]:
    print(g)

# Read GTF file
gtf_df = pd.read_csv(
    transcriptome_file,
    sep="\t",
    comment="#",
    header=None,
    names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
)

# Show example attribute column
print("\nFirst 5 raw attribute entries from GTF:")
print(gtf_df["attribute"].head())

# Extract gene_id
def extract_gene_id(attr):
    match = re.search(r'gene_id "([^"]+)"', attr)
    if match:
        return match.group(1)
    return None

gtf_df["gene_id"] = gtf_df["attribute"].apply(extract_gene_id).str.strip()

print("\nFirst 10 extracted gene IDs from GTF:")
print(gtf_df["gene_id"].dropna().head(10).tolist())

# --- Check intersection ---
matches = gene_ids.intersection(set(gtf_df["gene_id"]))
print(f"\nNumber of matching gene IDs: {len(matches)}")
print("Example matches:", list(matches)[:10])

# --- Filter ---
gtf_filtered = gtf_df[gtf_df["gene_id"].isin(gene_ids)][["gene_id", "seqname", "start", "end"]]
gtf_filtered.to_csv(output_file, sep="\t", index=False)
print(f"\nSaved {len(gtf_filtered)} gene locations to {output_file}")