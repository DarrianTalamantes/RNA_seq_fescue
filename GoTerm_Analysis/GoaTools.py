# This file will use GoaTools on the output of Deseq2 (2c_Deseq2_final_subset_VolcanoPlots.R and 2d_Upset_plots.R) 
# with the GoTerms idendified with Entap to make make a nice analysis

##########################################################
# Importing 
##########################################################

import os
import pandas as pd
import numpy as np
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudy

#Setting data directory
data = "/home/darrian/Documents/RNA_seq_fescue/Goatools_data"
feature_counts = "/home/darrian/Documents/RNA_seq_fescue/r_data/feature_counts.txt"
DEG_count_Data = "/home/darrian/Documents/RNA_seq_fescue/r_data/Endophytes_Up_Down_reg_HeatvsControl.csv" # Switch this line to change DEG set
Entap_Identified_Gos = "/home/darrian/Documents/RNA_seq_fescue/EnTAP_results/annotated_without_contam_gene_ontology_terms.tsv"
obo_loc=f"{data}/go-basic.obo"
os.makedirs(data, exist_ok=True)  

def main():
  
    
    # Importing the CSV file
    print("DIE DIE DIE")
    DEGs=import_csv(DEG_count_Data)
    DEGs.reset_index(inplace=True)
    DEGs.rename(columns={'index': 'Gene'}, inplace=True) 

    # Creating lists of just genes of different DEG groups
    # These are my study files
    filter_and_write_genes(DEGs,data_dir=data)

    # Getting list of all study files
    file_list = [f for f in os.listdir(data) if os.path.isfile(os.path.join(data, f))]

    # load in my assosiation file
    association_file = makeass_file(Entap_Identified_Gos, f"{data}/association.tsv")

    print("Assosiation file \n")
    for item in list(association_file)[:5]:
        print(item)

    #create or load in a population file 
    population_file = get_or_create_population_file(feature_counts, f"{data}/population.tsv")

    print("Population file \n")
    for item in list(population_file)[:5]:
        print(item)

    # Running goatools
    rungoatools(population_file, f"{data}/{file_list[1]}", association_file, obo_loc)
    print("loaded", f"{data}/{file_list[1]}")












def rungoatools(pop, study_loc, assoc, obo_loc):
    print("--------Running Goatools--------")

    # === Load Study ===
    with open(study_loc) as f:
        study = {line.strip() for line in f}

    # === Load GO DAG ===
    godag = GODag(obo_loc)

    goea = GOEnrichmentStudy(
        pop,    # All background genes
        assoc,         # Gene-to-GO associations
        godag,         # Ontologies
        methods=['fdr_bh']  # Multiple testing correction method
    )

    results = goea.run_study(study)

    print("\nSignificant GO terms (FDR < 0.05):")
    for r in results:
        if r.p_fdr_bh < 0.05:
            print(f"{r.GO}: {r.name} — p={r.p_fdr_bh:.4g} ({r.enrichment})")
            
    # goea.wr_tsv("goea_results.tsv", results)



def get_or_create_population_file(input_path, output_path):
    print("--------Population Function Running--------")
    """
    This creates the population file from annotated_without_contam Entap output
    If output_path exists, read and return it.
    Otherwise, create it by extracting the 'query_sequence' column from input_path,
    dropping duplicates, sorting, and writing to output_path.
    """
    if os.path.exists(output_path):
        # If the file already exists, just read it in
        print("Population file found, re-loading")
        with open(output_path) as f:
            df = {line.strip() for line in f}
        
    else:
        # Read input file
        df = pd.read_csv(input_path, sep="\t")
        df.rename(columns={df.columns[0]: "query_sequence"}, inplace=True)
        print(df.head())
        # Keep only the 'query_sequence' column
        df = df[["query_sequence"]]
        # Drop duplicates and sort
        df = df.drop_duplicates().sort_values(by="query_sequence")
        # Save to output_path without header or index
        df.to_csv(output_path, index=False, header=False)
        df = set(df["query_sequence"])
    return df



# This function looks for the assosiation file or makes it from a Entap output "annotated_without_contam"
def makeass_file(file_path, output_path="../Goatools_data/association.tsv"):
    print("--------Assosiation Function Running--------")
    # If the output file already exists, read it in.
    if os.path.exists(output_path):
        # File already exists, just read it (no headers)
        print("Assosiation file found, re-loading")
        assoc = {}
        with open(output_path) as f:
            for line in f:
                gene, go_terms = line.strip().split("\t", 1)
                assoc[gene] = set(go_terms.split(";"))
    else:
        # Create the association file
        df = pd.read_csv(file_path, sep="\t")
        df["go_id"] = df["go_id"].apply(lambda x: str(x).split("\t")[0])
        grouped = df.groupby("query_sequence")["go_id"] \
                    .apply(lambda x: ";".join(sorted(set(x)))) \
                    .reset_index()


        # Copy to avoid modifying original
        grouped = grouped.copy()

        # Trim gene IDs after 4th dot
        grouped.iloc[:, 0] = grouped.iloc[:, 0].apply(lambda x: ".".join(str(x).split(".")[:4]))

        # Remove rows that are duplicated (keep only fully unique rows)
        grouped = grouped[~grouped.duplicated(keep=False)]
        grouped.to_csv(output_path, sep="\t", index=False, header=False)


        # Now load the newly created file
        assoc = {}
        with open(output_path) as f:
            for line in f:
                gene, go_terms = line.strip().split("\t", 1)
                assoc[gene] = set(go_terms.split(";"))
        
    return assoc

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