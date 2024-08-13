# Purpose: This python script will be used to take lines from the tall fescue genome according to their gene_id
# Author: drt83172@uga.edu


def main():
    # Usage
    filter_gtf('/scratch/drt83172/Wallace_lab/RNA_SEQ/transcriptome/Fescue_transcriptome.gtf', 'significant_genes_list_unique.csv', '/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/significant_genes.gtf')

def filter_gtf(gtf_file, strings_file, output_file):
    # Read the list of strings from the file
    with open(strings_file, 'r') as f:
        strings = set(line.strip() for line in f)

    # Open the GTF file and the output file
    with open(gtf_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Check if the line contains any of the strings
            if any(s in line for s in strings):
                columns = line.split('\t')
                # Check if the third column is 'transcript'
                if len(columns) > 2 and columns[2] == 'transcript':
                    outfile.write(line)



if __name__ == '__main__':
    main()