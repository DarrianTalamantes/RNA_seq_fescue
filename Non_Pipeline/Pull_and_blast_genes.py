# Purpose: This python script will be used to take lines from the gtf file and seperate it into smaller gtf files. These then can be fed
# into the interproscan annotation pipline.
# Author: drt83172@uga.edu

import pandas as pd
import os
import glob


def main():
    # Take out all downregulated genes and signigficantly upregulated genes then put them into their ownfile
    significance_table = '/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/Non_Pipeline/significant_table.csv'  # Replace with the path to your input file
    dir_of_sig_files = '/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/significant_lists'   
    dir_of_gtfs = '/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/small_gtfs'
    split_downregulated_genes(significance_table, dir_of_sig_files)
    print("Splitting downregulated")
    split_upregulated_genes(significance_table, dir_of_sig_files)
    print("Splitting upregulated")


    # While loop to iterate over the file list
    file_list = list_files(dir_of_sig_files)
    index = 0
    while index < len(file_list):
        # Get the current file
        current_file = file_list[index]
        gtf_file = change_extension_to_gtf(current_file)
        current_list = dir_of_sig_files + '/' + current_file
        current_gtf_file = dir_of_gtfs + '/' + gtf_file
        filter_gtf('/scratch/drt83172/Wallace_lab/RNA_SEQ/transcriptome/Fescue_transcriptome.gtf', current_list, current_gtf_file)


# This function takes the significance table I made in the R file DeSeq2_Analysis.R and splits it into many smaller files.
def split_downregulated_genes(input_file, output_dir):
    # Read the data from the file
    df = pd.read_csv(input_file, sep="\t", index_col=0)  # Adjust separator if needed

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over each condition column (excluding the first column which is likely an index or identifier)
    for column in df.columns[1:]:
        # Filter for significantly downregulated genes
        downregulated = df[df[column].str.contains("Significant Downregulated", na=False)]
        
        # Define the output file name
        output_file = os.path.join(output_dir, f"{column}_significantly_downregulated.txt")
        
        # Save the filtered data to a new file
        downregulated[[column]].to_csv(output_file, sep="\t", index=True)

# This function takes the significance table I made in the R file DeSeq2_Analysis.R and splits it into many smaller files.
def split_upregulated_genes(input_file, output_dir):
    # Read the data from the file
    df = pd.read_csv(input_file, sep="\t", index_col=0)  # Adjust separator if needed

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over each condition column (excluding the first column which is likely an index or identifier)
    for column in df.columns[1:]:
        # Filter for significantly downregulated genes
        downregulated = df[df[column].str.contains("Significant Upregulated", na=False)]
        
        # Define the output file name
        output_file = os.path.join(output_dir, f"{column}_significantly_upregulated.txt")
        
        # Save the filtered data to a new file
        downregulated[[column]].to_csv(output_file, sep="\t", index=True)


def filter_gtf(gtf_file, strings_file, output_file):
    # Read the list of strings from the file
    with open(strings_file, 'r') as f:
        strings = set(line.strip() for line in f)

    # Dictionary to store lines by string
    string_lines = {s: [] for s in strings}

    # Open the GTF file and store lines that match
    with open(gtf_file, 'r') as infile:
        for line in infile:
            if any(s in line for s in strings):
                columns = line.split('\t')
                if len(columns) > 4 and columns[2] == 'transcript':
                    for s in strings:
                        if s in line:
                            string_lines[s].append(line)
                            break

    # Process each group of lines for each string
    with open(output_file, 'w') as outfile:
        for s, lines in string_lines.items():
            max_diff = None
            best_line = None
            
            for line in lines:
                columns = line.split('\t')
                start = int(columns[3])
                end = int(columns[4])
                diff = end - start

                if max_diff is None or diff > max_diff:
                    max_diff = diff
                    best_line = line

            if best_line:
                outfile.write(best_line)

# smal function to list everything in a directory
def list_files(directory):
    # List all files in the directory
    files = glob.glob(os.path.join(directory, '*'))
    # Extract just the filenames (optional)
    files = [os.path.basename(f) for f in files]
    return files

# Smol function to change extension to gtf
def change_extension_to_gtf(file_name):
    # Split the file name into name and extension
    base_name = os.path.splitext(file_name)[0]
    
    # Create the new file name with .gtf extension
    new_file_name = base_name + '.gtf'
    
    return new_file_name



if __name__ == '__main__':
    main()