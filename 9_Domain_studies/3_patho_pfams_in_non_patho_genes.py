### This script gets the pfam domains in the pathogenic genes and compares them with the pfam domains in the non-pathogenic genes
import pandas as pd

# Read the files
patho = pd.read_csv('./Domains/patho_unique_genes.txt', sep='\t')
print(patho.head())
non_patho = pd.read_csv('./Domains/non_patho_genes.txt', sep='\t')
print(non_patho.head())

# Get unique Pfam codes from pathogenic proteins
patho_pfam = set(patho['Pfam_found'].unique())

# Split the PFAM column in non_patho where multiple PFAMs are semicolon-separated
non_patho_expanded = non_patho['Pfam'].str.split(';', expand=True).stack().reset_index(drop=True)
non_patho_expanded = non_patho_expanded.str.strip()  # Remove any whitespace

# Get unique Pfam codes from pathogenic proteins
patho_pfam = set(patho['Pfam_found'].unique())

# Split the PFAM column in non_patho where multiple PFAMs are semicolon-separated
non_patho_expanded = non_patho['Pfam'].str.split(';', expand=True).stack().reset_index(drop=True)
non_patho_expanded = non_patho_expanded.str.strip()  # Remove any whitespace

# Find which Pfam codes from patho are present in non_patho
common_pfam = set(non_patho_expanded).intersection(patho_pfam)

# Count occurrences of each common Pfam in non_patho
pfam_counts = {}
for pfam in common_pfam:
    # Count in original patho dataframe
    patho_count = len(patho[patho['Pfam_found'] == pfam])
    # Count in non_patho (need to check original strings)
    non_patho_count = non_patho['Pfam'].str.contains(pfam, regex=False, na=False).sum()
    pfam_counts[pfam] = {'patho_count': patho_count, 'non_patho_count': non_patho_count}

# Create a result dataframe
result_df = pd.DataFrame.from_dict(pfam_counts, orient='index')
result_df.index.name = 'Pfam_code'
result_df = result_df.reset_index()

# Sort by Pfam code
result_df = result_df.sort_values('Pfam_code')

# Add a column for genes from both dataframes
patho_genes = patho.groupby('Pfam_found')['gene'].apply(lambda x: ';'.join(x)).to_dict()
non_patho_genes = non_patho.groupby('Pfam')['Gene Names'].apply(lambda x: ';'.join(x)).to_dict()

result_df['patho_genes'] = result_df['Pfam_code'].map(patho_genes)

# Modified function to handle NaN values
def get_non_patho_genes(pfam):
    # Fill NaN with empty string to avoid NA/NaN issues
    mask = non_patho['Pfam'].fillna('').str.contains(pfam, regex=False)
    matching_genes = non_patho[mask]['Gene Names']
    return ';'.join(matching_genes)

result_df['non_patho_genes'] = result_df['Pfam_code'].apply(get_non_patho_genes)

# Reorder columns
result_df = result_df[['Pfam_code', 'patho_count', 'non_patho_count', 'patho_genes', 'non_patho_genes']]

# Save to Excel
result_df.to_excel('./Domains/pfam_comparison.xlsx', index=False)

print("Excel file 'pfam_comparison.xlsx' has been created with the following columns:")
print("- Pfam_code: The PFAM identifier")
print("- patho_count: Number of occurrences in pathogenic proteins")
print("- non_patho_count: Number of occurrences in non-pathogenic proteins")
print("- patho_genes: Genes from pathogenic dataframe containing this Pfam")
print("- non_patho_genes: Genes from non-pathogenic dataframe containing this Pfam")