# Filter ClinVar processed file. We will conserve missense variants present in AD genes. Also, we will remove somatic variants.

import pandas as pd
import os


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


folder_path = "../2_Data_processing/data"
folder_path = os.path.normpath(folder_path)


# Set the path to the clinvar file
clinvar_file = f"{folder_path}/clinvar_homosapiens_processed.txt"
# Read data extracted from clinvar database using pandas package
clinvar = pd.read_csv(clinvar_file, delimiter="\t", header="infer")


folder_path = "../1_Data_collection/data"
folder_path = os.path.normpath(folder_path)
# Set the path to the file containing the AD genes from OMIM:
AD_OMIM_genes_file = f"{folder_path}/AD_genes_OMIM.txt"
# Open the file and retrieve the entries
AD_entries = open_file(AD_OMIM_genes_file)
print(type(AD_entries))
# Retrieve a list of AD genes
AD_genes = obtain_genes(AD_entries)

####################################################################################
##### IMPORTANT #####
### We noticed that are AD_genes when a deletion is produced but not for missense mutations.
# It means that we were conserving missense mutations as pathogenic for genes that are only AD when a deletion occurs.
# Mireia proposed a way to handle it: We conserve those AD genes that appear at least in three pathogenic variants
#####################

clinvar['Gene_name'] = clinvar['Gene(s)'].apply(lambda x: x.split('|')[0])

# First, we count the number of occurrences for each gene in the DataFrame
gene_counts = clinvar['Gene_name'].value_counts()

# Then, we create a new DataFrame that only includes genes that appear at least three times
genes_at_least_three = gene_counts[gene_counts >= 3]
print(genes_at_least_three)

# Now, we filter the original DataFrame to only include rows where the gene is in both AD_genes and genes_at_least_three
AD_clinvar = clinvar[clinvar["Gene(s)"].isin(AD_genes) & clinvar["Gene_name"].isin(genes_at_least_three.index)]
####################################################################################

# Save the entries in clinvar DataFrame that contain a gene in 'Gene(s)' column which is in the list of AD genes
#AD_clinvar = clinvar[clinvar["Gene(s)"].isin(AD_genes)]

# Filter empty rows for column Protein change
AD_clinvar = AD_clinvar[AD_clinvar["Protein change"].notna()]
AD_clinvar = AD_clinvar[AD_clinvar["Prot_change"].notna()]

AD_clinvar = AD_clinvar[AD_clinvar["Final aa"].notna()]
AD_clinvar = AD_clinvar[AD_clinvar["Initial aa"].notna()]


# Filter rows that contain a *
AD_clinvar = AD_clinvar[~AD_clinvar["Protein change"].str.contains("\*", regex=True)]

# Filter rows where the final aa or initial aa is: fs, Ter or del
exclude_cases = ["fs", "Ter", "del"]
AD_clinvar = AD_clinvar[~AD_clinvar["Final aa"].isin(exclude_cases)]
AD_clinvar = AD_clinvar[~AD_clinvar["Initial aa"].isin(exclude_cases)]

# Filter rows where the condition includes somatic cases
AD_clinvar = AD_clinvar[~AD_clinvar['Condition(s)'].str.contains('somatic', case=False)]

###
# Filter rows to avoid "Likely pathogenic" cases
#AD_clinvar = AD_clinvar[~AD_clinvar['Clinical significance (Last reviewed)'].str.contains('Likely pathogenic')]
###

# Set the output file name
output_file = "./data/clinvar_AD_missense_genes_min3.txt"
output_file = os.path.normpath(output_file)

# Save the new DataFrame in a csv file
AD_clinvar.to_csv(
    output_file,
    sep="\t",
    index=False,
)

# Print the results of the clinvar filtering step
print(
    f"From the {len(AD_entries)} AD genes, we kept all the names given by OMIM for the same gene, so we considered {len(AD_genes)} different gene names. The initial clinvar file had {len(clinvar)} entries. After filtering it for AD genes and missense variants, has {len(AD_clinvar)} entries. (We conserved AD genes that contained minimum 3 variants, to avoid those pathogenic AD genes but only for deletions.)"
)
