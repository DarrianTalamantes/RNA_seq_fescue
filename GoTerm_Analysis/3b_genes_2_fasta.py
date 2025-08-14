# Objective: This file will take the significant go term files and extract their gene-id and locations
# it will then pull the relevant fasta files so you can use BLAST on them.
# Results: This program finds that out of the results files that mattered we can not find anything in our transcriptome file that matches.
# Remember Scallop2 creates the transcriptome with the gene.XX.XX.XX names
# Transcriptome Final vs Transcriptome is the reason I cant find some gene.xx.xx.xx
 
import pandas as pd
import re
import unicodedata
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# --- Input file paths ---
go_terms_file = "../Goatools_data/Control_Down_results.txt"
transcriptome_file  = "/scratch/drt83172/Wallace_lab/RNA_SEQ/transcriptome/big/Fescue_transcriptome.gtf"
genome_file = "/scratch/drt83172/Wallace_lab/RNA_SEQ/Genome/Tall_fescue/tall_fescue_pv1.1.fasta"
output_file = "../Goatools_data/gene_locations.tsv"
output_fasta = "../Goatools_data/go2gene.fa"
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

# --- 1) Load GO terms & extract gene IDs ---
go_df = pd.read_csv(go_terms_file, sep="\t", dtype=str)
gene_ids = set()

if "study_items" in go_df.columns:
    for cell in go_df["study_items"].dropna():
        for g in re.findall(r'gene\.\d+(?:\.\d+)+', cell):
            gene_ids.add(norm(g))
else:
    with open(go_terms_file, "r", encoding="utf-8", errors="ignore") as fh:
        text = fh.read()
    for g in set(re.findall(r'gene\.\d+(?:\.\d+)+', text)):
        gene_ids.add(norm(g))

print(f"Loaded {len(gene_ids)} unique gene IDs from GO terms file")
print("Example gene IDs:")
for g in list(gene_ids)[:10]:
    print(" ", repr(g))

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