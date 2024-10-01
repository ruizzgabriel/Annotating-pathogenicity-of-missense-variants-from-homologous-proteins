# Filter Gnomad processed file. We will conserve missense variants present in AD genes.

import pandas as pd
import os
from subprocess import Popen, PIPE
import glob


def open_file(filename):
    """
    Open a txt file

    Args:
        filename(string): contains the name of the file with the path

    Returns:
        entries(list): containing all the entries of the file
    """
    # Open the file in reading mode
    with open(filename, "r") as file:
        # Retrieve the entries
        entries = file.readlines()
    # Close the file
    file.close()
    return entries


def obtain_genes(entries):
    """
    Obtain a list of genes from the entries of filtered OMIM file

    Args:
        entries(list): containing all the entries of the read file

    Returns:
        AD_genes(list): containing all genes in a more readable way
    """
    # Initialize the list of genes
    AD_genes = []
    # Iterate over the entries of the file
    for entry in entries:
        # Split the entries by commas
        entry = entry.split(",")
        # Iterate over the genes of each entry
        for gene in entry:
            # Strip the gene name
            gene = gene.strip()
            # Check if the gene is not present in the list of genes
            if gene not in AD_genes:
                # Add the name of the gene to the list of genes
                AD_genes.append(gene)
    return AD_genes


folder_path = "./data"
folder_path = os.path.normpath(folder_path)

# Set the path to the gnomAD file
#file_pattern = "gnomad_r*_processed.csv"
#file_pattern = "gnomad_r*_processed_rescuedErrors.csv"
file_pattern = "gnomad_r4_GRCh38_AD_missense_genes_joined.txt"    ############# !!!!!!!!!!!!!
# Search the files matching the pattern
gnomad_file_list = glob.glob(f"{folder_path}/{file_pattern}")
print(gnomad_file_list)

folder_path = "../2_Data_processing/data"
folder_path = os.path.normpath(folder_path)
# Set the path to the file containing the AD genes from OMIM:
AD_OMIM_genes_file = f"{folder_path}/AD_genes_OMIM.txt"
# Open the file and retrieve the entries
AD_entries = open_file(AD_OMIM_genes_file)
# Retrieve a list of AD genes
AD_genes = obtain_genes(AD_entries)

# Iterate over the different version files of gnomAD
for gnomad_file in gnomad_file_list:
    if "r2" in gnomad_file:
        gnomad_v = "gnomad_r2_1"
        genome = "GRCh37"
    elif "r3" in gnomad_file:
        gnomad_v = "gnomad_r3"
        genome = "GRCh38"
    elif "r4" in gnomad_file:
        gnomad_v = "gnomad_r4"
        genome = "GRCh38"

    # Read the first line (header)
    header = pd.read_csv(gnomad_file, nrows=1)
    gnomad_df = pd.read_csv(gnomad_file, header=None, skiprows=1)  # Ignore first line
    gnomad_df.columns = header.columns  # Add the header

    # Filter by 'missense_variant'
    gnomad_df = gnomad_df[gnomad_df['consequence'] == 'missense_variant']
    #print(gnomad_df)
    
    # Save the entries in gnomAD DataFrame that contain a gene in 'gene' column which is in the list of AD genes
    AD_gnomad = gnomad_df[gnomad_df["gene"].isin(AD_genes)]

    # Filter entries by genome and exome allele frequency                    ############# !!!!!!!!!!!!!
    AD_gnomad = AD_gnomad[(AD_gnomad['genome_af'] > 10e-6) | (AD_gnomad['exome_af'] > 10e-6)]  # e.g. 10e-5,  5*10e-4

    # There are rows in gnomAD in which we are not interested --> we remove them
    values = ['Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic', 'Pathogenic; drug response', 'Pathogenic; other']
    AD_gnomad = AD_gnomad[~AD_gnomad['clinical_significance'].isin(values)]

    # Select rows where 'clinical_significance' do not contain 'Uncertain significance' and similar cases.
    AD_gnomad = AD_gnomad[~AD_gnomad['clinical_significance'].astype(str).str.contains('Uncertain significance')]
    AD_gnomad = AD_gnomad[~AD_gnomad['clinical_significance'].astype(str).str.contains('Conflicting interpretations of pathogenicity')]
    AD_gnomad = AD_gnomad[~AD_gnomad['clinical_significance'].astype(str).str.contains('Likely pathogenic')]

    # Set the output file name
    output_file = f"./data/{gnomad_v}_{genome}_AD_missense_genes_freq-6_joined.txt"          ############# !!!!!!!!!!!!!
    output_file = os.path.normpath(output_file)

    # Save the new DataFrame in a csv file
    AD_gnomad.to_csv(output_file, sep=",", index=False)
