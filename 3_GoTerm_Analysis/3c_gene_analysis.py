# THis file will take the gene file(s) made from 2d_Upset_plots.R and extract their sequences
# The idea here is to find the function of these genes. 





import pandas as pd                   # this must be installed with ml not conda now
import re
import unicodedata
from Bio import SeqIO                 # this must be installed with ml not conda now
from Bio.SeqRecord import SeqRecord

# --- Input file paths ---
gene_list = "../Gene_lists/shared_genes_heat_heatxpercipitation.txt"
transcriptome_file  = "/scratch/drt83172/Wallace_lab/RNA_SEQ/transcriptome/big/Fescue_transcriptome.gtf"
genome_file = "/scratch/drt83172/Wallace_lab/RNA_SEQ/Genome/Tall_fescue/tall_fescue_pv1.1.fasta"
output_file = "../Gene_lists/shared_heatxheatpercipitation.tsv"
output_fasta = "../Gene_lists/shared_heatxheatpercipitation.fa"
debug_gene = "gene.41791.92.0"   # set to any gene you want to sanity check

# --- Helpers ---
def norm(s):
    """Normalize strings to avoid hidden differences (BOM, zero-width, CRLF, extra spaces)."""
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return None
    s = str(s)
    s = unicodedata.normalize("NFKC", s)
    s = s.replace("\ufeff", "").replace("\u200b", "")  # BOM, zero-width space
    return s.strip()

# --- 1) Load gene IDs ---
with open(gene_list) as f:
    gene_ids = set(line.strip() for line in f if line.strip())

print(f"Loaded {len(gene_ids)} gene IDs")
print(list(gene_ids)[:10])  # show first 10

# --- 2) Load GTF & extract gene_id ---
gtf_df = pd.read_csv(
    transcriptome_file,
    sep="\t",
    comment="#",
    header=None,
    names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"],
    dtype={"seqname": str, "source": str, "feature": str, "start": int, "end": int,
           "score": str, "strand": str, "frame": str, "attribute": str}
)

def extract_gene_id(attr):
    if attr is None:
        return None
    m = re.search(r'gene_id "([^"]+)"', attr)
    return m.group(1) if m else None

gtf_df["gene_id"] = gtf_df["attribute"].apply(extract_gene_id).apply(norm)

print("\nExample 'attribute' → extracted gene_id from GTF:")
for _, row in gtf_df.head(5).iterrows():
    print(f"ATTR: {repr(row['attribute'])}  --> gene_id: {repr(row['gene_id'])}")

# --- 3) Intersection check ---
go_set = set(filter(None, (norm(g) for g in gene_ids)))
gtf_set = set(filter(None, gtf_df["gene_id"].tolist()))
matches = go_set.intersection(gtf_set)

print(f"\nIntersection size: {len(matches)} (GO ∩ GTF)")
print("Example matches:", [repr(x) for x in list(matches)[:10]])

if debug_gene:
    print(f"\nDEBUG check for {debug_gene!r}:")
    print("  in GO set?  ", debug_gene in go_set)
    print("  in GTF set? ", debug_gene in gtf_set)
    if debug_gene in gtf_set:
        print("  GTF rows for this gene:")
        print(gtf_df.loc[gtf_df["gene_id"] == debug_gene, ["seqname","start","end"]].head(3))

# --- 4) Filter matches ---
gtf_filtered = gtf_df.loc[gtf_df["gene_id"].isin(go_set), ["gene_id", "seqname", "start", "end", "strand"]].drop_duplicates()

print("\nExample matches found in transcriptome:")
for _, r in gtf_filtered.head(10).iterrows():
    print(f"{r['gene_id']} -> {r['seqname']}:{r['start']}-{r['end']} ({r['strand']})")

gtf_filtered.to_csv(output_file, sep="\t", index=False)
print(f"\nSaved {len(gtf_filtered)} gene locations to {output_file}")

# --- 5) Extract sequences from genome ---
print("\nLoading genome into memory...")
genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

records = []
for _, row in gtf_filtered.iterrows():
    seqname = row["seqname"]
    start = row["start"]
    end = row["end"]
    strand = row["strand"]
    gene_id = row["gene_id"]

    if seqname not in genome_dict:
        print(f"WARNING: {seqname} not found in genome, skipping {gene_id}")
        continue

    # GTF is 1-based inclusive; Biopython slices are 0-based exclusive at end
    seq = genome_dict[seqname].seq[start-1:end]
    if strand == "-":
        seq = seq.reverse_complement()

    records.append(SeqRecord(seq, id=gene_id, description=f"{seqname}:{start}-{end}({strand})"))

SeqIO.write(records, output_fasta, "fasta")
print(f"Saved {len(records)} sequences to {output_fasta}")