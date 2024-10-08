# Purpose: This python script will be used to take lines from the gtf file and seperate it into smaller gtf files. These then can be fed
# into the interproscan annotation pipline.
# Author: drt83172@uga.edu

import pandas as pd
import os
import glob
import re

def main():
    # Take out all downregulated genes and signigficantly upregulated genes then put them into their ownfile
    significance_table = '/scratch/drt83172/Wallace_lab/RNA_SEQ/Scripts/RNA_seq_fescue/Non_Pipeline/significant_table.csv'  # Replace with the path to your input file
    dir_of_sig_files = '/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/significant_lists'   
    dir_of_gtfs = '/scratch/drt83172/Wallace_lab/RNA_SEQ/manual_annotation/small_gtfs'
    
    '''
    split_downregulated_genes(significance_table, dir_of_sig_files)
    print("Splitting downregulated")
    split_upregulated_genes(significance_table, dir_of_sig_files)
    print("Splitting upregulated")

    
    # While loop to iterate over the file list
    # This while loop subsets the big gtf file to smaller ones based on the significance lists
    file_list = list_files(dir_of_sig_files)

    index = 0
    while index < len(file_list):
        # Get the current file
        current_file = file_list[index]
        gtf_file = change_extension_to_gtf(current_file)
        current_list = dir_of_sig_files + '/' + current_file
        current_gtf_file = dir_of_gtfs + '/' + gtf_file
        print("making the file ", gtf_file)
        filter_gtf('/scratch/drt83172/Wallace_lab/RNA_SEQ/transcriptome/Fescue_transcriptome.gtf', current_list, current_gtf_file)
        index += 1
    '''

    '''
    # This function iterates through all the small gtf files and removes duplicate genes from them.
    file_list = list_files(dir_of_gtfs)
    index = 0
    while index < len(file_list):
        current_file = file_list[index]
        print("processing " + current_file)
        dir_and_current_file = dir_of_gtfs + '/' + current_file
        current_file2 = dir_of_gtfs + "/dupped_" + current_file
        gtf_dup_remover(dir_and_current_file,current_file2 )
        index += 1
    '''


# searches for genes in big gtf file and subsets to make a small one of chosen genes from strings_file
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

# This takes a gtf file as an input and then removes the duplicated gene ids.
def gtf_dup_remover(input_file, output_file):
    # Read the data from the file
    df = pd.read_csv(input_file, sep="\t", header=None, engine='python')

    # Extract the gene_id from the last column
    df['gene_id'] = df[8].apply(lambda x: re.search(r'gene_id "([^"]+)"', x).group(1) if re.search(r'gene_id "([^"]+)"', x) else None)

    # Calculate the difference between column 5 and column 4
    df['difference'] = df[4] - df[3]

    # Drop duplicates based on gene_id, keeping the row with the maximum difference
    result = df.loc[df.groupby('gene_id')['difference'].idxmax()]

    # Drop unnecessary columns and keep the required ones
    result = result[[0, 1, 2, 3, 4, 5, 6, 'gene_id', 'difference']]

    # Save the processed DataFrame to a new file
    result.to_csv(output_file, sep="\t", header=False, index=False)


# This function takes the significance table I made in the R file DeSeq2_Analysis.R and splits it into many smaller files.
def split_downregulated_genes(input_file, output_dir):
    # Read the data from the file
    df = pd.read_csv(input_file, sep=",", index_col=0)  # Adjust separator if needed

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over each condition column (excluding the first column which is likely an index or identifier)
    for column in df.columns[2:]:
        # Filter for significantly downregulated genes
        downregulated = df[df[column].str.contains("Significant Downregulated", na=False)]
        
        # Define the output file name
        output_file = os.path.join(output_dir, f"{column}_significantly_downregulated.txt")
        
        # Only keep the first column (index or identifier) and drop the second column
        downregulated_genes = downregulated.index

        # Save the filtered data to a new file
        pd.DataFrame(downregulated_genes, columns=['Gene']).to_csv(output_file, sep="\t", index=False)

# This function takes the significance table I made in the R file DeSeq2_Analysis.R and splits it into many smaller files.
def split_upregulated_genes(input_file, output_dir):
    # Read the data from the file
    df = pd.read_csv(input_file, sep=",", index_col=0)  # Adjust separator if needed

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over each condition column (excluding the first column which is likely an index or identifier)
    for column in df.columns[2:]:
        # Filter for significantly downregulated genes
        downregulated = df[df[column].str.contains("Significant Upregulated", na=False)]
        
        # Define the output file name
        output_file = os.path.join(output_dir, f"{column}_significantly_upregulated.txt")
       
        # Only keep the first column (index or identifier) and drop the second column
        downregulated_genes = downregulated.index

        # Save the filtered data to a new file
        pd.DataFrame(downregulated_genes, columns=['Gene']).to_csv(output_file, sep="\t", index=False)

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