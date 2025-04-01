### This script reads all the downloaded file from uniprot which contains all homo sapiens
### genes and it's associated pfam codes and selects the non-pathogenic ones. For selecting
### the non-pathogenic ones we will use the list of pathogenic genes and the others that are not pathogneic
### (at least 3 missense pathogenic mutatons in clinvar) will be selcted

import pandas as pd

patho = pd.read_csv('./12_Domain_studies/Domains/patho_unique_genes.txt', sep = '\t')
#print(patho.head())
print(len(patho))

non_patho = pd.read_csv('./12_Domain_studies/Domains/20250401_homo_sapiens_genes_and_pfams.tsv', sep = '\t')
#print(non_patho.head())

# Drop rows where 'Gene Names' or 'Pfam' is NaN
non_patho = non_patho.dropna(subset=['Gene Names', 'Pfam'])

# Split 'Gene Names' by space and select the first element
non_patho['Gene Names'] = non_patho['Gene Names'].str.split().str[0]
print(len(non_patho))
#print(non_patho.head())

# Filter non_patho to keep only rows where Gene Names is not in patho's gene column
non_patho_filtered = non_patho[~non_patho['Gene Names'].isin(patho['gene'])]
print(len(non_patho_filtered))

# Drop the 'Entry' column
non_patho_filtered = non_patho_filtered.drop('Entry', axis=1)

# Save to text file with tab separation
non_patho_filtered.to_csv('./12_Domain_studies/Domains/non_patho_filtered.txt', sep='\t', index=False)

print("File saved as 'non_patho_filtered.txt'")
print(non_patho_filtered)