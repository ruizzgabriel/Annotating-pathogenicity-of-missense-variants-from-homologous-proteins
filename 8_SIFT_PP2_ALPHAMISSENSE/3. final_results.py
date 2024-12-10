import pandas as pd
import numpy as np


def add_verify_sift_column(df):
    # Define conditions
    conditions = [
        (df['SIFT_prediction'] == 'Deleterious') & (df['db'] == 'gnomAD'),
        (df['SIFT_prediction'] == 'Neutral') & (df['db'] == 'ClinVar'),
        (df['SIFT_prediction'] == 'Deleterious') & (df['db'] == 'ClinVar'),
        (df['SIFT_prediction'] == 'Neutral') & (df['db'] == 'gnomAD'),
        (df['SIFT_prediction'] == ' ')
    ]

    # Define corresponding values
    values = [0, 0, 1, 1, np.nan]

    # Apply conditions to create the 'Verify_SIFT' column
    df['Verify_SIFT'] = np.select(conditions, values, default=np.nan)
    return df

def add_verify_polyphen_column(df):
    # Define conditions
    conditions = [
        (df['POLY_prediction'] == 'probably damaging') & (df['db'] == 'gnomAD'),
        (df['POLY_prediction'] == 'possibly damaging') & (df['db'] == 'gnomAD'),
        (df['POLY_prediction'] == 'benign') & (df['db'] == 'gnomAD'),
        (df['POLY_prediction'] == 'probably damaging') & (df['db'] == 'ClinVar'),
        (df['POLY_prediction'] == 'possibly damaging') & (df['db'] == 'ClinVar'),
        (df['POLY_prediction'] == 'benign') & (df['db'] == 'ClinVar'),
        (df['POLY_prediction'] == 'unknown'),
        (df['POLY_prediction'] == ' ')
    ]

    # Define corresponding values
    values = [0, 0, 1, 1, 1, 0, np.nan, np.nan]

    # Apply conditions to create the 'Verify_SIFT' column
    df['Verify_POLYPHEN'] = np.select(conditions, values, default=np.nan)
    return df




# Create a mapping dictionary from 1-letter to 3-letter amino acid codes
amino_acid_mapping = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr'
}


uniprot = pd.read_csv('uniprotkb_AND_reviewed_true_AND_model_o_2024_10_28.tsv', sep = '\t')

#################################
###### STRICT PAIRS  ############
#################################
df_sp = pd.read_excel('Vars_from_strict_pairs.xlsx')

df_sp_sift = pd.read_csv('./SIFT_results/Strict_pairs_SIFT_result.tsv', sep = '\t')
df_sp_sift = df_sp_sift[['PREDICTION (cutoff=-2.5)']]
# Rename a single column
df_sp_sift.rename(columns={'PREDICTION (cutoff=-2.5)': 'SIFT_prediction'}, inplace=True)


df_sp_poly = pd.read_csv('./POLYPHEN2_results/Strict_pairs_POLYPHEN2_result.txt', sep = '\t')
df_sp_poly = df_sp_poly[['#o_acc               ', ' o_pos', 'o_aa1', 'o_aa2', '        prediction']]
df_sp_poly.rename(columns={'#o_acc               ' : 'accession', ' o_pos' : 'position', 'o_aa1' : 'aa1', 'o_aa2' : 'aa2', '        prediction' : 'POLY_prediction'}, inplace=True)
#print(df_sp, df_sp_sift, df_sp_poly)



# Merge the two DataFrames based on the condition
df_sp = df_sp.merge(uniprot[['Entry Name', 'Entry']], how='left', left_on='uniprot', right_on='Entry Name')
# Optionally, you can drop the 'Entry Name' column if it's no longer needed
df_sp = df_sp.drop(columns=['Entry Name'])



#### SIFT
# Join SIFT result
df_sp = df_sp.join(df_sp_sift, how='left')
# Verify SIFT result
df_sp = add_verify_sift_column(df_sp)


#### POLYPHEN2
# Trim spaces and newline characters from all cells in df_sp_poly
df_sp_poly = df_sp_poly.applymap(lambda x: x.strip() if isinstance(x, str) else x)

# Convert 1-letter codes to 3-letter codes in df_sp_poly
df_sp_poly['aa1'] = df_sp_poly['aa1'].replace(amino_acid_mapping)
df_sp_poly['aa2'] = df_sp_poly['aa2'].replace(amino_acid_mapping)

# Convert the 'position' column to integers
df_sp['pos_prot_change'] = df_sp['pos_prot_change'].astype(float)

df_sp_poly = df_sp_poly.drop_duplicates(subset=['accession', 'aa1', 'position', 'aa2'])

# Perform a left join
df_sp = df_sp.merge(
    df_sp_poly[['accession', 'position', 'aa1', 'aa2', 'POLY_prediction']],  # Select only the necessary columns from df_sp_poly
    left_on=['Entry', 'initial_aa', 'pos_prot_change', 'final_aa'],  # Columns from df_sp
    right_on=['accession', 'aa1', 'position', 'aa2'],  # Columns from df_sp_poly
    how='left'  # Left join
)

df_sp.drop(columns=['accession', 'aa1', 'position', 'aa2'], inplace=True)

# Verify POLYPHEN result
df_sp = add_verify_polyphen_column(df_sp)


df_sp.to_csv('./Final_results/Evaluated_vars_from_strict_pairs.txt', index=False)

#print(df_sp)


###########################################
###### SIMILAR FINAL AA PAIRS  ############
###########################################

df_sf = pd.read_excel('Vars_from_similar_final_AA_pairs.xlsx')

df_sf_sift = pd.read_csv('./SIFT_results/Similar_final_aa_pairs_SIFT_result.tsv', sep = '\t')
df_sf_sift = df_sf_sift[['PREDICTION (cutoff=-2.5)']]
# Rename column
df_sf_sift.rename(columns={'PREDICTION (cutoff=-2.5)': 'SIFT_prediction'}, inplace=True)


df_sf_poly = pd.read_csv('./POLYPHEN2_results/Similar_final_aa_pairs_POLYPHEN2_result.txt', sep = '\t')
df_sf_poly = df_sf_poly[['#o_acc               ', ' o_pos', 'o_aa1', 'o_aa2', '        prediction']]
df_sf_poly.rename(columns={'#o_acc               ' : 'accession', ' o_pos' : 'position', 'o_aa1' : 'aa1', 'o_aa2' : 'aa2', '        prediction' : 'POLY_prediction'}, inplace=True)
#print(df_sp, df_sp_sift, df_sf_poly)



# Merge the two DataFrames based on the condition
df_sf = df_sf.merge(uniprot[['Entry Name', 'Entry']], how='left', left_on='uniprot', right_on='Entry Name')
# Drop the 'Entry Name'
df_sf = df_sf.drop(columns=['Entry Name'])



#### SIFT
# Join SIFT result
df_sf = df_sf.join(df_sf_sift, how='left')
# Verify SIFT result
df_sf = add_verify_sift_column(df_sf)


#### POLYPHEN2
# Trim spaces and newline characters from all cells in df_sf_poly
df_sf_poly = df_sf_poly.applymap(lambda x: x.strip() if isinstance(x, str) else x)

# Convert 1-letter codes to 3-letter codes in df_sf_poly
df_sf_poly['aa1'] = df_sf_poly['aa1'].replace(amino_acid_mapping)
df_sf_poly['aa2'] = df_sf_poly['aa2'].replace(amino_acid_mapping)

# Convert the 'position' column to integers
df_sf['pos_prot_change'] = df_sf['pos_prot_change'].astype(float)

df_sf_poly = df_sf_poly.drop_duplicates(subset=['accession', 'aa1', 'position', 'aa2'])

# Perform a left join
df_sf = df_sf.merge(
    df_sf_poly[['accession', 'position', 'aa1', 'aa2', 'POLY_prediction']],  # Select only the necessary columns from df_sp_poly
    left_on=['Entry', 'initial_aa', 'pos_prot_change', 'final_aa'],  # Columns from df_sp
    right_on=['accession', 'aa1', 'position', 'aa2'],  # Columns from df_sp_poly
    how='left'  # Left join
)

df_sf.drop(columns=['accession', 'aa1', 'position', 'aa2'], inplace=True)

# Verify POLYPHEN result
df_sf = add_verify_polyphen_column(df_sf)


df_sf.to_csv('./Final_results/Evaluated_vars_from_similar_final_AA_pairs.txt', index=False)




#############################################
###### SIMILAR INITIAL AA PAIRS  ############
#############################################

df_si = pd.read_excel('Vars_from_similar_initial_AA_pairs.xlsx')

df_si_sift = pd.read_csv('./SIFT_results/Similar_initial_aa_pairs_SIFT_result.tsv', sep = '\t')
df_si_sift = df_si_sift[['PREDICTION (cutoff=-2.5)']]
# Rename column
df_si_sift.rename(columns={'PREDICTION (cutoff=-2.5)': 'SIFT_prediction'}, inplace=True)


df_si_poly = pd.read_csv('./POLYPHEN2_results/Similar_initial_aa_pairs_POLYPHEN2_result.txt', sep = '\t')
df_si_poly = df_si_poly[['#o_acc               ', ' o_pos', 'o_aa1', 'o_aa2', '        prediction']]
df_si_poly.rename(columns={'#o_acc               ' : 'accession', ' o_pos' : 'position', 'o_aa1' : 'aa1', 'o_aa2' : 'aa2', '        prediction' : 'POLY_prediction'}, inplace=True)
#print(df_sp, df_sp_sift, df_sf_poly)



# Merge the two DataFrames based on the condition
df_si = df_si.merge(uniprot[['Entry Name', 'Entry']], how='left', left_on='uniprot', right_on='Entry Name')
# Drop the 'Entry Name' column
df_si = df_si.drop(columns=['Entry Name'])



#### SIFT
# Join SIFT result
df_si = df_si.join(df_si_sift, how='left')
# Verify SIFT result
df_si = add_verify_sift_column(df_si)


#### POLYPHEN2
# Trim spaces and newline characters from all cells in df_sf_poly
df_si_poly = df_si_poly.applymap(lambda x: x.strip() if isinstance(x, str) else x)

# Convert 1-letter codes to 3-letter codes in df_sf_poly
df_si_poly['aa1'] = df_si_poly['aa1'].replace(amino_acid_mapping)
df_si_poly['aa2'] = df_si_poly['aa2'].replace(amino_acid_mapping)

# Convert the 'position' column to integers
df_si['pos_prot_change'] = df_si['pos_prot_change'].astype(float)

df_si_poly = df_si_poly.drop_duplicates(subset=['accession', 'aa1', 'position', 'aa2'])

# Perform a left join
df_si = df_si.merge(
    df_si_poly[['accession', 'position', 'aa1', 'aa2', 'POLY_prediction']],  # Select only the necessary columns from df_sp_poly
    left_on=['Entry', 'initial_aa', 'pos_prot_change', 'final_aa'],  # Columns from df_sp
    right_on=['accession', 'aa1', 'position', 'aa2'],  # Columns from df_sp_poly
    how='left'  # Left join
)

df_si.drop(columns=['accession', 'aa1', 'position', 'aa2'], inplace=True)

# Verify POLYPHEN result
df_si = add_verify_polyphen_column(df_si)


df_si.to_csv('./Final_results/Evaluated_vars_from_similar_initial_AA_pairs.txt', index=False)




#####################################################
###### SIMILAR INITIAL + FINAL AA PAIRS  ############
#####################################################

df_sif = pd.read_excel('Vars_from_similar_i+f_AA_pairs.xlsx')

df_sif_sift = pd.read_csv('./SIFT_results/Similar_i+f_aas_pairs_SIFT_result.tsv', sep = '\t')
df_sif_sift = df_sif_sift[['PREDICTION (cutoff=-2.5)']]
# Rename column
df_sif_sift.rename(columns={'PREDICTION (cutoff=-2.5)': 'SIFT_prediction'}, inplace=True)


df_sif_poly = pd.read_csv('./POLYPHEN2_results/Similar_i+f_aas_pairs_POLYPHEN2_result.txt', sep = '\t')
df_sif_poly = df_sif_poly[['#o_acc               ', ' o_pos', 'o_aa1', 'o_aa2', '        prediction']]
df_sif_poly.rename(columns={'#o_acc               ' : 'accession', ' o_pos' : 'position', 'o_aa1' : 'aa1', 'o_aa2' : 'aa2', '        prediction' : 'POLY_prediction'}, inplace=True)
#print(df_sp, df_sp_sift, df_sf_poly)



# Merge the two DataFrames based on the condition
df_sif = df_sif.merge(uniprot[['Entry Name', 'Entry']], how='left', left_on='uniprot', right_on='Entry Name')
# Drop the 'Entry Name' column 
df_sif = df_sif.drop(columns=['Entry Name'])



#### SIFT
# Join SIFT result
df_sif = df_sif.join(df_sif_sift, how='left')
# Verify SIFT result
df_sif = add_verify_sift_column(df_sif)


#### POLYPHEN2
# Trim spaces and newline characters from all cells in df_sf_poly
df_sif_poly = df_sif_poly.applymap(lambda x: x.strip() if isinstance(x, str) else x)

# Convert 1-letter codes to 3-letter codes in df_sf_poly
df_sif_poly['aa1'] = df_sif_poly['aa1'].replace(amino_acid_mapping)
df_sif_poly['aa2'] = df_sif_poly['aa2'].replace(amino_acid_mapping)

# Convert the 'position' column to integers
df_sif['pos_prot_change'] = df_sif['pos_prot_change'].astype(float)

df_sif_poly = df_sif_poly.drop_duplicates(subset=['accession', 'aa1', 'position', 'aa2'])

# Perform a left join
df_sif = df_sif.merge(
    df_sif_poly[['accession', 'position', 'aa1', 'aa2', 'POLY_prediction']],  # Select only the necessary columns from df_sp_poly
    left_on=['Entry', 'initial_aa', 'pos_prot_change', 'final_aa'],  # Columns from df_sp
    right_on=['accession', 'aa1', 'position', 'aa2'],  # Columns from df_sp_poly
    how='left'  # Left join
)

df_sif.drop(columns=['accession', 'aa1', 'position', 'aa2'], inplace=True)

# Verify POLYPHEN result
df_sif = add_verify_polyphen_column(df_sif)


df_sif.to_csv('./Final_results/Evaluated_vars_from_similar_i+f_AA_pairs.txt', index=False)