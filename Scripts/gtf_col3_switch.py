# This script takes a gtf file and put the geneID in column three so bamtools getfasta will put it in the name of its output

import sys
import re

def swap_gene_id(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.strip():
                fields = line.strip().split('\t')
                match = re.search(r'gene_id "([^"]+)"', fields[8])
                if match:
                    gene_id = match.group(1)
                    fields[2] = gene_id
                    outfile.write('\t'.join(fields) + '\n')

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    swap_gene_id(input_file, output_file)