"""
This script processes UniProt data for Homo sapiens genes and their associated Pfam codes.
It separates the genes into pathogenic and non-pathogenic categories based on a provided
list of pathogenic proteins, and saves both sets to separate files.
"""

import pandas as pd

# File paths
PATHO_FILE = './12_Domain_studies/Domains/Pathogenic AD proteins.xlsx'
UNIPROT_FILE = './12_Domain_studies/Domains/20250401_homo_sapiens_genes_and_pfams.tsv'
PATHO_OUTPUT = './12_Domain_studies/Domains/patho_genes.txt'
NON_PATHO_OUTPUT = './12_Domain_studies/Domains/non_patho.txt'

def load_and_preview_pathogenic_genes(file_path):
    """Load and display initial information about pathogenic genes."""
    df = pd.read_excel(file_path)
    print("Pathogenic Genes Preview:")
    print(df.head())
    print(f"Total number of pathogenic genes: {len(df)}")
    return df

def load_and_clean_uniprot(file_path):
    """Load UniProt data and clean it by removing NaN values and splitting gene names."""
    df = pd.read_csv(file_path, sep='\t')
    # Remove rows with missing 'Gene Names' or 'Pfam' values
    df = df.dropna(subset=['Gene Names', 'Pfam'])
    # Take the first gene name if multiple are listed
    df['Gene Names'] = df['Gene Names'].str.split().str[0]
    print(f"Total genes after cleaning: {len(df)}")
    return df

def split_pathogenic_nonpathogenic(uniprot_df, patho_df):
    """Split UniProt data into pathogenic and non-pathogenic based on pathogenic gene list."""
    # Filter for pathogenic genes
    patho_genes = uniprot_df[uniprot_df['Gene Names'].isin(patho_df['Pathogenic_proteins'])]
    # Filter for non-pathogenic genes
    non_patho = uniprot_df[~uniprot_df['Gene Names'].isin(patho_df['Pathogenic_proteins'])]
    
    print(f"Number of pathogenic genes found: {len(patho_genes)}")
    print(f"Number of non-pathogenic genes: {len(non_patho)}")
    return patho_genes, non_patho

def save_results(patho_df, non_patho_df, patho_output, non_patho_output):
    """Save pathogenic and non-pathogenic dataframes to files after dropping 'Entry' column."""
    # Remove 'Entry' column from both dataframes
    patho_df = patho_df.drop('Entry', axis=1)
    non_patho_df = non_patho_df.drop('Entry', axis=1)
    
    # Save to tab-separated files
    patho_df.to_csv(patho_output, sep='\t', index=False)
    non_patho_df.to_csv(non_patho_output, sep='\t', index=False)
    
    print(f"Files saved as '{patho_output}' and '{non_patho_output}'")

def main():
    """Main function to execute the gene separation process."""
    # Load data
    patho = load_and_preview_pathogenic_genes(PATHO_FILE)
    uniprot = load_and_clean_uniprot(UNIPROT_FILE)
    
    # Split into pathogenic and non-pathogenic
    patho_genes, non_patho = split_pathogenic_nonpathogenic(uniprot, patho)
    
    # Save results
    save_results(patho_genes, non_patho, PATHO_OUTPUT, NON_PATHO_OUTPUT)
    
    # Display samples
    print("\nPathogenic Genes Sample:")
    print(patho_genes.head())
    print("\nNon-pathogenic Genes Sample:")
    print(non_patho.head())

if __name__ == "__main__":
    main()