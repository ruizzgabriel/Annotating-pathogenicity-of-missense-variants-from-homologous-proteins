#### This script gets 2 files:
### File 1: patho_unique_genes.txt containing the unique genes with Pfam domains from pathogenic variants
### File 2: non_patho_unique_genes.txt containing the unique genes from non-pathogenic variants

import pandas as pd

# Read the Excel file
excel_file = pd.ExcelFile('./Domains/Supplementary_material.xlsx')

# Convert tabs to DataFrames
patho = excel_file.parse('S4.Vars_mapped_SEED_Pfam_aligns')
non_patho = excel_file.parse('S3.Non-patho missense variants')

# Get unique genes with Pfam for pathogenic variants
# We'll keep both Gene(s) and Pfam columns and remove duplicates
patho_unique_genes = patho[['gene', 'Pfam_found']].drop_duplicates()

# Get unique genes for non-pathogenic variants (unchanged)
non_patho_unique_genes = pd.DataFrame(
    non_patho['gene'].unique(),
    columns=['Gene(s)']
)

# Save each DataFrame to a separate text file
patho_unique_genes.to_csv('./Domains/patho_unique_genes.txt', index=False, sep='\t')
non_patho_unique_genes.to_csv('./Domains/non_patho_unique_genes.txt', index=False, sep='\t')

# Print confirmation
print("Saved unique genes with Pfam domains from pathogenic variants to 'patho_unique_genes.txt'")
print(f"Number of unique gene-Pfam combinations: {len(patho_unique_genes)}")
print("Saved unique genes from non-pathogenic variants to 'non_patho_unique_genes.txt'")
print(f"Number of unique genes: {len(non_patho_unique_genes)}")