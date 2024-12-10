import pandas as pd

df_entries = pd.read_csv('uniprotkb_AND_reviewed_true_AND_model_o_2024_10_28.tsv', sep = '\t')
#print(df_entries)

def add_accession_column(df, df_entries):
    # Create a dictionary from df_entries to map 'Entry Name' to 'Entry'
    entry_dict = df_entries.set_index('Entry Name')['Entry'].to_dict()
    
    # Use the dictionary to map 'uniprot' values in df_sp to 'Entry' values
    df['accession'] = df['uniprot'].map(entry_dict)
    
    return df


def select_columns(df, name):
    # Select the desired columns
    name = df[['accession', 'pos_prot_change', 'initial_aa', 'final_aa']]
    return name

aa_mapping = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
    'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
    'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
    'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
    'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
}

def convert_aa(df):
    # Convert 'initial_aa' and 'final_aa' columns
    df['initial_aa'] = df['initial_aa'].map(aa_mapping)
    df['final_aa'] = df['final_aa'].map(aa_mapping)
    return df



#################################
###### STRICT PAIRS  ############
#################################
df_sp = pd.read_excel('Vars_from_strict_pairs.xlsx')
#print(df_sp)

df_sp = add_accession_column(df_sp, df_entries)
df_sp = select_columns(df_sp, df_sp)
df_sp = convert_aa(df_sp)

df_sp.to_csv('Query_vars_from_strict_pairs.txt', sep=" ", index=False, header=False)

###########################################
###### SIMILAR FINAL AA PAIRS  ############
###########################################
df_sf = pd.read_excel('Vars_from_similar_final_AA_pairs.xlsx')
#print(df_sf)

df_sf = add_accession_column(df_sf, df_entries)
df_sf = select_columns(df_sf, df_sf)
df_sf = convert_aa(df_sf)

df_sf.to_csv('Query_vars_from_similar_final_AA_pairs.txt', sep=" ", index=False, header=False)



#############################################
###### SIMILAR INITIAL AA PAIRS  ############
#############################################
df_si = pd.read_excel('Vars_from_similar_initial_AA_pairs.xlsx')
#print(df_si)

df_si = add_accession_column(df_si, df_entries)
df_si = select_columns(df_si, df_si)
df_si = convert_aa(df_si)

df_si.to_csv('Query_vars_from_similar_initial_AA_pairs.txt', sep=" ", index=False, header=False)



#####################################################
###### SIMILAR INITIAL + FINAL AA PAIRS  ############
#####################################################
df_sif = pd.read_excel('Vars_from_similar_i+f_AA_pairs.xlsx')
#print(df_sif)

df_sif = add_accession_column(df_sif, df_entries)
df_sif = select_columns(df_sif, df_sif)
df_sif = convert_aa(df_sif)

df_sif.to_csv('Query_vars_from_similar_i+f_AA_pairs.txt', sep=" ", index=False, header=False)




print("Go to the following pages and use to_polyphen_sift.txt")
print("- Polyphen2: http://genetics.bwh.harvard.edu/pph2/bgi.shtml")
print(
    "- SIFT: http://sift.jcvi.org/protein_batch_submit.php?species=human PROVEAN Protein Batch"
)
print("From Polyphen2 download the SHORT results file!")
