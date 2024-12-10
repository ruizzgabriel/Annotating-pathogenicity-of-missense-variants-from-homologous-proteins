import pandas as pd

# Dictionary to map one-letter amino acid codes to three-letter codes
aa_dict = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    'B': 'Asx', 'Z': 'Glx', 'X': 'Xaa', '*': 'Ter'
}

# Function to convert one-letter to three-letter codes
def convert_to_three_letter(aa):
    return aa_dict.get(aa, 'Unknown')



#################################
###### STRICT PAIRS  ############
#################################

df =  pd.read_excel('Strict_pairs.xlsx')
#print(df)

# Create two separate dataframes, one for each set of values
df1 = df[['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 'gene_name', 'pos_prot_change', 'db1', 'uniprot1', 'score_substmat', 'seq_identity(%)', 'Pair_type']]
df2 = df[['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 'gene_name2', 'pos_prot_change2', 'db2', 'uniprot2', 'score_substmat', 'seq_identity(%)', 'Pair_type']]

# Rename columns in df2 to match the structure of df1
df2.columns = ['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 'gene_name', 'pos_prot_change', 'db1', 'uniprot1', 'score_substmat', 'seq_identity(%)', 'Pair_type']

# Combine the two dataframes by alternating their rows
result = pd.concat([df1, df2], axis=0).sort_index(kind='merge').reset_index(drop=True)

# Remove rows where 'gene_name' or 'pos_prot_change' is missing
result = result.dropna(subset=['gene_name', 'pos_prot_change'])

# Rename columns in result 
result.columns = ['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 'gene_name', 'pos_prot_change', 'db', 'uniprot', 'score_substmat', 'seq_identity(%)', 'Pair_type']

result.to_excel('Vars_from_strict_pairs.xlsx', index=False)



###########################################
###### SIMILAR FINAL AA PAIRS  ############
###########################################

df =  pd.read_excel('Similar_final_aa_pairs.xlsx')
#print(df)

# Create two separate dataframes, one for each set of values
df1 = df[['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 'gene_name', 'pos_prot_change', 'db1', 'uniprot1', 'score_substmat', 'seq_identity(%)', 'Pair_type']]
df2 = df[['Pfam_code', 'initial_aa', 'pos_align', 'final_aa2', 'gene_name2', 'pos_prot_change2', 'db2', 'uniprot2', 'score_substmat', 'seq_identity(%)', 'Pair_type']]

# Rename columns in df2 to match the structure of df1
df2.columns = ['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 'gene_name', 'pos_prot_change', 'db1', 'uniprot1', 'score_substmat', 'seq_identity(%)', 'Pair_type']

# Combine the two dataframes by alternating their rows
result = pd.concat([df1, df2], axis=0).sort_index(kind='merge').reset_index(drop=True)

# Remove rows where 'gene_name' or 'pos_prot_change' is missing
result = result.dropna(subset=['gene_name', 'pos_prot_change'])

# Rename columns in result 
result.columns = ['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 'gene_name', 'pos_prot_change', 'db', 'uniprot', 'score_substmat', 'seq_identity(%)', 'Pair_type']

# Apply the conversion to the 'final_aa' column
result['final_aa'] = result['final_aa'].apply(convert_to_three_letter)

result.to_excel('Vars_from_similar_final_AA_pairs.xlsx', index=False)



#############################################
###### SIMILAR INITIAL AA PAIRS  ############
#############################################

df =  pd.read_excel('Similar_initial_aa_pairs.xlsx')
#print(df)

# Create two separate dataframes, one for each set of values
df1 = df[['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 'gene_name', 'pos_prot_change', 'db1', 'uniprot1', 'score_substmat', 'seq_identity(%)', 'Pair_type']]
df2 = df[['Pfam_code', 'initial_aa2', 'pos_align', 'final_aa', 'gene_name2', 'pos_prot_change2', 'db2', 'uniprot2', 'score_substmat', 'seq_identity(%)', 'Pair_type']]

# Rename columns in df2 to match the structure of df1
df2.columns = ['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 'gene_name', 'pos_prot_change', 'db1', 'uniprot1', 'score_substmat', 'seq_identity(%)', 'Pair_type']

# Combine the two dataframes by alternating their rows
result = pd.concat([df1, df2], axis=0).sort_index(kind='merge').reset_index(drop=True)

# Remove rows where 'gene_name' or 'pos_prot_change' is missing
result = result.dropna(subset=['gene_name', 'pos_prot_change'])

# Rename columns in result 
result.columns = ['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 'gene_name', 'pos_prot_change', 'db', 'uniprot', 'score_substmat', 'seq_identity(%)', 'Pair_type']

# Apply the conversion to the 'initial_aa' column
result['initial_aa'] = result['initial_aa'].apply(convert_to_three_letter)

result.to_excel('Vars_from_similar_initial_AA_pairs.xlsx', index=False)



#####################################################
###### SIMILAR INITIAL + FINAL AA PAIRS  ############
#####################################################

df =  pd.read_excel('Similar_i+f_aas_pairs.xlsx')
#print(df)

# Create two separate dataframes, one for each set of values
df1 = df[['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 'gene_name', 'pos_prot_change', 'db1', 'uniprot1', 'score_substmat', 'seq_identity(%)', 'Pair_type']]
df2 = df[['Pfam_code', 'initial_aa2', 'pos_align', 'final_aa2', 'gene_name2', 'pos_prot_change2', 'db2', 'uniprot2', 'score_substmat', 'seq_identity(%)', 'Pair_type']]

# Rename columns in df2 to match the structure of df1
df2.columns = ['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 'gene_name', 'pos_prot_change', 'db1', 'uniprot1', 'score_substmat', 'seq_identity(%)', 'Pair_type']

# Combine the two dataframes by alternating their rows
result = pd.concat([df1, df2], axis=0).sort_index(kind='merge').reset_index(drop=True)

# Remove rows where 'gene_name' or 'pos_prot_change' is missing
result = result.dropna(subset=['gene_name', 'pos_prot_change'])

# Rename columns in result 
result.columns = ['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 'gene_name', 'pos_prot_change', 'db', 'uniprot', 'score_substmat', 'seq_identity(%)', 'Pair_type']

# Apply the conversion to the 'initial_aa' column
result['initial_aa'] = result['initial_aa'].apply(convert_to_three_letter)
result['final_aa'] = result['final_aa'].apply(convert_to_three_letter)

result.to_excel('Vars_from_similar_i+f_AA_pairs.xlsx', index=False)