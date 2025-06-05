# This file will use GoaTools on the output of Deseq2 (2c_Deseq2_final_subset_VolcanoPlots.R and 2d_Upset_plots.R) 
# with the GoTerms idendified with Entap to make make a nice analysis

##########################################################
# Importing 
##########################################################

import os
import pandas as pd
import numpy as np
import goatools

#Setting data directory
data = "/home/darrian/Documents/RNA_seq_fescue/tempdata"
os.makedirs(data, exist_ok=True)  

def main():
  
    
    # Importing the CSV file
    DEG_count_Data="/home/darrian/Documents/RNA_seq_fescue/r_data/Treatments_Up_Down_reg.csv"
    print("DIE DIE DIE")
    DEGs=import_csv(DEG_count_Data)
    DEGs.reset_index(inplace=True)
    DEGs.rename(columns={'index': 'Gene'}, inplace=True) 

    # Displaying the data and fixing it if necessary
    print(DEGs.head())
    filter_and_write_genes(DEGs,data_dir=data)



    # Running goatools
    rungoatools(1,2,3,4)

def rungoatools(pop, study, assoc, go_dag):
    print("This Ran LOL")



def filter_and_write_genes(df, gene_col="Gene", output_prefix="genes", data_dir= data):
    """
    For each column in the DataFrame (except the gene column), filter out rows
    where the column is '0', and write the gene names of remaining rows to a file.
    
    Parameters:
    - df: Pandas DataFrame
    - gene_col: name of the column containing gene names (default: "gene")
    - output_prefix: prefix for output files (default: "genes")
    """
    for col in df.columns:
        if col == gene_col:
            continue  # Skip the gene column
        
        # Filter out rows where the current column is '0' (as string or number)
        filtered = df[df[col] != "0"]
        filtered = filtered[filtered[col] != 0]  # Also check numeric zero
        
        # Extract gene names
        gene_list = filtered[gene_col]
        
        # Save to file
        output_filename = f"{output_prefix}.{col}.txt"
        gene_list.to_csv(f"{data_dir}/{output_filename}", index=False, header=False)
        print(f"Wrote {len(gene_list)} genes to {output_filename}")


def import_csv(filepath):
    """
    Imports a CSV file and returns a pandas DataFrame.

    Parameters:
    filepath (str): Path to the CSV file.

    Returns:
    pd.DataFrame: DataFrame containing the CSV data.
    """
    try:
        df = pd.read_csv(filepath,index_col=0)
        print(f"Successfully loaded file: {filepath}")
        return df
    except FileNotFoundError:
        print(f"File not found: {filepath}")
    except pd.errors.EmptyDataError:
        print("File is empty.")
    except pd.errors.ParserError:
        print("Error parsing the file.")
main()