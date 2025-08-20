# This will take a fasta file and BLAST every sequence in the file. Good to check any genes we want to know the functions of.


#!/usr/bin/env python3
import os
import subprocess
from Bio import SeqIO

# Input FASTA
fasta_file = "/home/darrian/Documents/RNA_seq_fescue/Gene_lists/shared_heatxheatpercipitation.fa"
output_file  = "/home/darrian/Documents/RNA_seq_fescue/Gene_lists/shared_heatxheatpercipitation_BLASTED.txt"
filtered_fasta_path = "/home/darrian/Documents/RNA_seq_fescue/Gene_lists/shared_heatxheatpercipitation_longest.fa"


def main():
    filtered_fasta = filter_longest_per_gene(fasta_file, filtered_fasta_path)
    
    # Open one output file for all top hits
    with open(output_file, "w", buffering=1) as out_handle:
        # Loop through each sequence
        for record in SeqIO.parse(filtered_fasta, "fasta"):
            seq_id = record.id
            tmp_seq_file = f"{seq_id}.tmp.fasta"
            
            # Write sequence to temporary FASTA
            with open(tmp_seq_file, "w") as f:
                f.write(f">{seq_id}\n{str(record.seq)}\n")
            
            # Run remote BLAST against NCBI nt with descriptive output
            cmd = [
                "blastn",
                "-query", tmp_seq_file,
                "-db", "nt",
                "-remote",
                "-outfmt", "6 qseqid sseqid pident length evalue bitscore stitle",
                "-evalue", "1e-5",
                "-max_target_seqs", "1"  # Only get top hit
            ]
            
            print(f"Running remote BLAST for {seq_id}...")
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Write results with divider
            out_handle.write(f"\n{'='*80}\n")
            out_handle.write(f"Top hit for {seq_id}\n")
            out_handle.write(f"{'='*80}\n")
            out_handle.write(result.stdout.strip() + "\n")
            
            # Cleanup temporary FASTA
            os.remove(tmp_seq_file)

    print(f"\n✅ All BLAST searches completed! Top hits are in: {output_file}")



def filter_longest_per_gene(input_fasta, output_fasta):
    """
    Keep only the longest sequence for each gene in a FASTA file.
    
    Parameters:
    - input_fasta: str, path to the original FASTA
    - output_fasta: str, path to save the filtered FASTA
    """
    longest_seqs = {}

    for record in SeqIO.parse(input_fasta, "fasta"):
        # Extract geneID (everything before first space)
        gene_id = record.id.split()[0]

        # Keep the longest sequence per gene
        if gene_id not in longest_seqs or len(record.seq) > len(longest_seqs[gene_id].seq):
            longest_seqs[gene_id] = record

    # Write filtered sequences to output
    with open(output_fasta, "w") as out_fh:
        SeqIO.write(longest_seqs.values(), out_fh, "fasta")

    print(f"Filtered FASTA saved to {output_fasta} ({len(longest_seqs)} sequences kept)")
    return output_fasta

main()
    