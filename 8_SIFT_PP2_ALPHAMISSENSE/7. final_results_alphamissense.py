import pandas as pd
import os

# Load the datasets
sp = pd.read_csv('./Final_results/Evaluated_vars_from_strict_pairs.txt')
si = pd.read_csv('./Final_results/Evaluated_vars_from_similar_initial_AA_pairs.txt')
sf = pd.read_csv('./Final_results/Evaluated_vars_from_similar_final_AA_pairs.txt')
sif = pd.read_csv('./Final_results/Evaluated_vars_from_similar_i+f_AA_pairs.txt')

# Mapping from 3-letter to 1-letter amino acid codes
aa_3to1 = {
    'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E', 'Phe': 'F',
    'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Lys': 'K', 'Leu': 'L',
    'Met': 'M', 'Asn': 'N', 'Pro': 'P', 'Gln': 'Q', 'Arg': 'R',
    'Ser': 'S', 'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y'
}

def convert_aa(aa_3):
    """Convert three-letter amino acid code to one-letter code."""
    return aa_3to1.get(aa_3, '')

def add_alphamissense_prediction_general(df, folder_path, 
                                         initial_aa_col, pos_col, final_aa_col, entry,
                                         prediction_col='ALPHAMISSENSE_prediction'):
    """
    Adds ALPHAMISSENSE prediction to any DataFrame with specified column names.
    
    Parameters:
    - df (pd.DataFrame): The input DataFrame.
    - folder_path (str): Path to the folder containing ENST files.
    - initial_aa_col (str): Name of the column with initial amino acid (3-letter code).
    - pos_col (str): Name of the column with position (can be float).
    - final_aa_col (str): Name of the column with final amino acid (3-letter code).
    - entry (str): Name of the column with the ENST transcript ID.
    - prediction_col (str): Name of the new column to store the ALPHAMISSENSE prediction (default: 'ALPHAMISSENSE_prediction').
    
    Returns:
    - pd.DataFrame: Updated DataFrame with the ALPHAMISSENSE prediction added.
    """
    # Initialize the new prediction column with None values
    df[prediction_col] = None
    
    # Iterate through each row in df
    for index, row in df.iterrows():
        entry_id = row[entry]
        initial_aa_3 = row[initial_aa_col]
        pos = int(row[pos_col])
        final_aa_3 = row[final_aa_col]
        
        # Convert 3-letter amino acids to 1-letter codes
        initial_aa_1 = convert_aa(initial_aa_3)
        final_aa_1 = convert_aa(final_aa_3)
        
        # Construct the mutation string to search (e.g., "M1A")
        mutation = f"{initial_aa_1}{pos}{final_aa_1}"
        
        # Define the file path for the ENST file
        file_path = os.path.join(folder_path, f"{entry_id}.txt")
        
        # Check if the file exists
        if os.path.exists(file_path):
            # Read the file and search for the mutation
            with open(file_path, 'r') as file:
                for line in file:
                    parts = line.strip().split("\t")
                    if len(parts) >= 4 and parts[1] == mutation:
                        # If mutation is found, add the pathogenicity to the prediction column
                        df.at[index, prediction_col] = parts[3]
                        break
    return df

# Folder path 
folder_path = "./Alphamissense/Alphamissense_by_accessions"

# Dictionary of dataframes and file paths for easy iteration
datasets = {
    "sp": ("./Final_results/Evaluated_vars_from_strict_pairs(alphamissense).txt", sp),
    "si": ("./Final_results/Evaluated_vars_from_similar_initial_AA_pairs(alphamissense).txt", si),
    "sf": ("./Final_results/Evaluated_vars_from_similar_final_AA_pairs(alphamissense).txt", sf),
    "sif": ("./Final_results/Evaluated_vars_from_similar_i+f_AA_pairs(alphamissense).txt", sif)
}

# Apply the function to each dataframe and save the result
for name, (file_path, df) in datasets.items():
    updated_df = add_alphamissense_prediction_general(
        df, folder_path,
        initial_aa_col='initial_aa', 
        pos_col='pos_prot_change', 
        final_aa_col='final_aa',
        entry='Entry'
    )
    # Save the updated dataframe to the specified file path
    updated_df.to_csv(file_path, index=False)
    print(f"Processed and saved {name} dataset to {file_path}.")
