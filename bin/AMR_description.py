#!/usr/bin/env python3
# Description: This script add custom AMR description file to abricate output file. 
#The script takes three arguments: the path to the AMR custom description file, the path to the input abricate ouput CSV file, and the name of the output file. 

import sys

amr_description_file = sys.argv[1]  # Path to the AMR description file
csv_file = sys.argv[2]  # Path to the input CSV file
outfile_name = sys.argv[3]  # Name of the output file

custom_AMR_gene = {}  # Dictionary to store custom AMR gene descriptions

# Read the AMR description file and populate the custom_AMR_gene dictionary
for line in open(amr_description_file, "r"):
    line = line.strip()
    gene, description = line.split("\t")
    custom_AMR_gene[gene] = description

output_file = f"{outfile_name}_CARD.txt"  # Path to the output file

# Open the output file for writing
with open(output_file, "w") as f:
    # Read the input CSV file line by line
    for line in open(csv_file, "r"):
        line = line.strip()
        if line.startswith("#"):
            # Print and write the header line as it is
            f.write(line + "\n")
        else:
            parts = line.split("\t")
            gene = parts[5]  # Extract the gene name from the CSV line
            if gene in custom_AMR_gene:
                # If the gene is present in the custom_AMR_gene dictionary, update the description
                parts[13] = custom_AMR_gene[gene]
            # Print and write the modified CSV line
                f.write("\t".join(parts) + "\n")




        
           