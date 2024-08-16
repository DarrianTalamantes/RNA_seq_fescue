# Purpose: This script will look at all the outputs of interpro scan and create a csv file that counts how many times each protein appears.

import os
import csv
import pandas as pd
from collections import defaultdict

# Directory containing the TSV files
input_dir = "/home/darrian/Desktop/UGA/Wallace_Lab/RNA_seq_fescue/interpro_results"

# Initialize a dictionary to store counts of each unique value
unique_value_counts = defaultdict(lambda: defaultdict(int))

# Get a list of all TSV files in the directory
tsv_files = [f for f in os.listdir(input_dir) if f.endswith('.tsv')]

# Process each TSV file
for tsv_file in tsv_files:
    file_path = os.path.join(input_dir, tsv_file)
    
    # Read the TSV file
    df = pd.read_csv(file_path, sep='\t', header=None)
    
    # Extract column 6 (index 5) and count occurrences of each value
    for value in df[5].dropna():
        unique_value_counts[value][tsv_file] += 1

# Create a list of all unique values
all_unique_values = sorted(unique_value_counts.keys())

# Prepare the header for the CSV file (first column + file names)
header = ['Unique_Value'] + tsv_files

# Write the results to a CSV file
output_file = "proteinCount.csv"
with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)
    
    for value in all_unique_values:
        row = [value] + [unique_value_counts[value][tsv_file] for tsv_file in tsv_files]
        writer.writerow(row)

print(f"Output written to {output_file}")